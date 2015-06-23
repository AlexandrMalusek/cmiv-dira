/**
 *  Arrays are indexed differently when copying from MATLAB code.
 *  This means that an 2x3 array will be indexed:
 *  0 2 4
 *  1 3 5
 *
 *  instead of row-wise first:
 *  0 1 2
 *  3 4 5
 *
 */

#include <math.h>
#include "mex.h"

static void 
MD2(double *Wei2Ptr, double *densPtr, double *atte1matPtr, double *atte2matPtr,
    double *att2Ptr, double *dens2Ptr, bool *maskPtr, int M,
    int N, int image_size);

static char rcs_id[] = "$Revision: 1.10 $";

/* Input Arguments */
#define ATTE1MAT  (prhs[0])
#define ATTE2MAT  (prhs[1])
#define ATT2      (prhs[2])
#define DENS2     (prhs[3])
#define MASK      (prhs[4])

/* Output Arguments */
#define	WEI2      (plhs[0])
#define DENS      (plhs[1])

/**
 * Need 5 input arguments: AttE1mat, AttE2mat, Att2, Dens2, mask
 * Produces 2 outputs: Wei2, dens
 */

void 
mexFunction(int nlhs, mxArray  *plhs[], int nrhs, const mxArray  *prhs[])
{
  double *atte1matPtr;  /* 511x511 pointer to array of measured attenuation coefficients at energy E1 */
  double *atte2matPtr;  /* 511x511 pointer to array of measured attenuation coefficients at energy E2 */
  double *att2Ptr;    /* 2x2 pointer to array of tabled attenuation coefficients*/
  double *dens2Ptr;     /* 1x2 pointer to mass density of doublet base materials */
  bool *maskPtr;      /* 511x511 pointer to mask defining the tissue to be decomposed*/
  
  /* temp pointers for copying values */
  double *from_doublePtr;   
  double *to_doublePtr;   
  bool *boolPtr1;
  bool *boolPtr2;
  
  int N;    /* columns of the image to process */
  int M;     /* rows of the image to process */
  int image_size;     /* Total image size, columns * rows */
  
  int att2_columns;
  int att2_rows;
  int att2_size;
  
  int dens2_columns;
  int dens2_rows;
  int dens2_size;
  
  int k;          /* Loop counter */
  int ndim;         /* Number of dimensions for WEI2 output array */
  int dims[3];      /* Array specifying the dimensions for WEI2 array */
  
  /* Check validity of arguments */
  if (nrhs != 5)
  {
    mexErrMsgTxt("Incorrect number of INPUT arguments.");
  }
  
  if (nlhs != 2)
  {
    mexErrMsgTxt("Incorrect number of OUTPUT arguments.");
  }
  
  if (mxIsSparse(ATTE1MAT) || mxIsSparse(ATTE2MAT) || mxIsSparse(ATT2) || mxIsSparse(DENS2) || mxIsSparse(MASK))
  {
    mexErrMsgTxt("Sparse inputs not supported.");
  }
  
  if (!mxIsDouble(ATTE1MAT) || !mxIsDouble(ATTE2MAT) || !mxIsDouble(ATT2) || !mxIsDouble(DENS2) )
  {
    mexErrMsgTxt("Input must be double.");
  }
  
  if (!mxIsLogical(MASK))
  {
    mexErrMsgTxt("MASK must be logical.");
  }
  
  /* Allocate memory for the arrays */
 
  N = mxGetN(ATTE1MAT);
  M = mxGetM(ATTE1MAT);
  image_size = M* N; 
  atte1matPtr = (double *) mxCalloc(image_size, sizeof(double));
  atte2matPtr = (double *) mxCalloc(image_size, sizeof(double));
  maskPtr = (bool *) mxCalloc(image_size, sizeof(bool));
  
  att2_columns = mxGetN(ATT2);
  att2_rows = mxGetM(ATT2);
  att2_size = att2_rows * att2_columns;
  att2Ptr = (double *) mxCalloc(att2_size, sizeof(double));
  
  dens2_columns = mxGetN(DENS2);
  dens2_rows = mxGetM(DENS2);
  dens2_size = dens2_rows* dens2_columns; 
  dens2Ptr = (double *) mxCalloc(dens2_size, sizeof(double));
  
  
  /* Copy all data */
  from_doublePtr = mxGetPr(ATTE1MAT);
  to_doublePtr = atte1matPtr;
  for (k = 0; k < image_size; k++)
  *(to_doublePtr++) = *(from_doublePtr++);
  
  from_doublePtr = mxGetPr(ATTE2MAT);
  to_doublePtr = atte2matPtr;
  for (k = 0; k < image_size; k++)
  *(to_doublePtr++) = *(from_doublePtr++);
  
  from_doublePtr = mxGetPr(ATT2);
  to_doublePtr = att2Ptr;
  for (k = 0; k < att2_size; k++)
  *(to_doublePtr++) = *(from_doublePtr++);
  
  from_doublePtr = mxGetPr(DENS2);
  to_doublePtr = dens2Ptr;
  for (k = 0; k < dens2_size; k++)
  *(to_doublePtr++) = *(from_doublePtr++);
  
  boolPtr1 = mxGetLogicals(MASK);
  boolPtr2 = maskPtr;
  for (k = 0; k < image_size; k++)
  *(boolPtr2++) = *(boolPtr1++);
  
  DENS = mxCreateDoubleMatrix(M, N, mxREAL);
  
  ndim = 3;
  dims[0] = M;
  dims[1] = N;
  dims[2] = 2;
  WEI2 = mxCreateNumericArray(ndim, dims, mxDOUBLE_CLASS, mxREAL);
  
  MD2(mxGetPr(WEI2), mxGetPr(DENS), atte1matPtr, atte2matPtr, att2Ptr, dens2Ptr,
      maskPtr,M, N, image_size);
}

static void 
MD2(double *Wei2Ptr, double *densPtr, double *atte1matPtr, double *atte2matPtr,
    double *att2Ptr, double *dens2Ptr, bool *maskPtr, int M,
    int N, int image_size)
{  
  /* Loop variables */
  int x,y;
  
  /* Matrices for linear equation */
  double m[2][2];
  double b[2];
  double w[2];
  double quota;  
  
  for(y=0;y<M;++y)
  {    
    for(x=0;x<N;++x)
    {
      if(maskPtr[y*M + x] == 1)
      {
        m[0][0] = att2Ptr[0]/dens2Ptr[0] - att2Ptr[2]/dens2Ptr[1];
        m[0][1] = -atte1matPtr[y*M + x];
        m[1][0] = att2Ptr[1]/dens2Ptr[0] - att2Ptr[3]/dens2Ptr[1];
        m[1][1] = -atte2matPtr[y*M + x];
        
        b[0] = -att2Ptr[2]/dens2Ptr[1];
        b[1] = -att2Ptr[3]/dens2Ptr[1];
        
        w[0] = 0;
        w[1] = 0;
        
        /* Use Gaussian elimination to solve the linear equation */
        quota   = m[1][0]/m[0][0];
        m[1][1] = m[1][1] - (m[0][1]*quota);
        b[1]    = b[1] - (b[0]*quota);
        w[1]    = b[1]/(m[1][1]);
        b[0]    = b[0]-(m[0][1]*w[1]);
        w[0]    = b[0]/(m[0][0]);

        Wei2Ptr[y*M + x] = w[0];
        Wei2Ptr[y*M + x + image_size] = 1 - w[0];         
        densPtr[y*M + x] = 1/w[1]; 
      }
    }
  } 
}   