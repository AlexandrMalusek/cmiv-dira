/*-------------------------------------------------------------------*/
/* This routine computes a Joseph sinogram from a pixelized phantom. */
/* In a Joseph sinogram, the interpolation kernel is designed        */
/* relative to the phantom grid.                                     */
/* Different interpolation filters are available.                    */
/* The routine is based on Matlabs iradon.m.                         */
/*                                                                   */
/* Written by Maria Magnusson Seger 2003-04                          */
/*-------------------------------------------------------------------*/

#include <math.h>
#include "mex.h"

static void 
MD2(double *Wei2Ptr, double *densPtr, double *atte1matPtr, double *atte2matPtr,
    double *att2Ptr, double *dens2Ptr, bool *maskPtr, int image_rows,
    int image_columns, int image_size);

static char rcs_id[] = "$Revision: 1.10 $";

#define MAXX(x,y) ((x) > (y) ? (x) : (y))

#define Pi 3.14159265358979
#define SMALL 0.000000001

/* Input Arguments */
#define ATTE1MAT  (prhs[0])
#define ATTE2MAT  (prhs[1])
#define ATT2      (prhs[2])
#define DENS2     (prhs[3])
#define MASK      (prhs[4])

/* Output Arguments */
#define WEI2      (plhs[0])
#define DENS      (plhs[1])

/**
 *  nlhs  Number of output (left-side) arguments, or the size of the plhs array.
 *  plhs  Array of output arguments.
 *  nrhs  Number of input (right-side) arguments, or the size of the prhs array.
 *  prhs  Array of input arguments.
 *
 */


/**
 * Need 5 input arguments: AttE1mat, AttE2mat, Att2, Dens2, mask
 * Produces 2 outputs: Wei2, dens
 */

void 
mexFunction(int nlhs, mxArray  *plhs[], int nrhs, const mxArray  *prhs[])
{
  double *atte1matPtr;  /* 511x511 pointer to array of measured attenuation coefficients at energy E1 */
  double *atte2matPtr;  /* 511x511 pointer to array of measured attenuation coefficients at energy E2 */
  double *att2Ptr;      /* 2x2 pointer to array of tabled attenuation coefficients*/
  double *dens2Ptr;     /* 1x2 pointer to mass density of doublet base materials */
  bool *maskPtr;	/* 511x511 pointer to mask defining the tissue to be decomposed*/
  
  /* temp pointers for copying values */
  double *from_doublePtr;   
  double *to_doublePtr;   
  bool *boolPtr1;
  bool *boolPtr2;
  
  int image_columns;  /* columns of the image to process */
  int image_rows;     /* rows of the image to process */
  int image_size;     /* Total image size, columns * rows */
  
  int att2_columns;
  int att2_rows;
  int att2_size;
  
  int dens2_columns;
  int dens2_rows;
  int dens2_size;
  
  int k;            /* Loop counter */
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
 
  image_columns = mxGetN(ATTE1MAT);
  image_rows = mxGetM(ATTE1MAT);
  image_size = image_rows* image_columns; 
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
  
  DENS = mxCreateDoubleMatrix(image_rows, image_columns, mxREAL);
  
  ndim = 3;
  dims[0] = image_rows;
  dims[1] = image_columns;
  dims[2] = 3;
  WEI2 = mxCreateNumericArray(ndim, dims, mxDOUBLE_CLASS, mxREAL);
  
  MD2(mxGetPr(WEI2), mxGetPr(DENS), atte1matPtr, atte2matPtr, att2Ptr, dens2Ptr,
      maskPtr,image_rows, image_columns, image_size);
}

static void 
MD2(double *Wei2Ptr, double *densPtr, double *atte1matPtr, double *atte2matPtr,
    double *att2Ptr, double *dens2Ptr, bool *maskPtr, int image_rows,
    int image_columns, int image_size)
{    
  /* Loop variables */
  int k;
  int l;
  
  /* Matrices for linear equation */
  double M[2][2];
  double b[2];
  double w[2];
  double quota;    
  
  for(k = 0; k < image_rows; ++k)
    {        
      for(l = 0; l < image_columns; ++l)
        {
	  if(*(maskPtr+(k*image_rows)+l) == 1)
            {
	      M[0][0] = (*att2Ptr) / (*dens2Ptr) - (*(att2Ptr+2)) / (*(dens2Ptr+1));
	      M[0][1] = -(*(atte1matPtr+(k*image_rows)+l));
	      M[1][0] = (*(att2Ptr+1) / (*dens2Ptr)) - (*(att2Ptr+3)) / (*(dens2Ptr+1));
	      M[1][1] = -(*(atte2matPtr+(k*image_rows)+l));
              
	      b[0] = -(*(att2Ptr+2)) / (*(dens2Ptr+1));
	      b[1] = -(*(att2Ptr+3)) / (*(dens2Ptr+1));
              
	      w[0] = 0;
	      w[1] = 0;
              
	      /* Use Gaussian elimination to solve the linear equation */
	      quota = M[1][0] / M[0][0];
	      M[1][1] = M[1][1] - (M[0][1] * quota);
	      b[1] = b[1] - (b[0] * quota);
	      w[1] = b[1] / (M[1][1]);
	      b[0] = b[0] - (M[0][1] * w[1]);
	      w[0] = b[0] / (M[0][0]);
	      
	      *(Wei2Ptr + k*image_rows + l) = w[0];
	      *(Wei2Ptr + k*image_rows + l + image_size) = 1.0 - w[0];               
	      *(densPtr + k*image_rows + l) = 1.0 / w[1]; 
            }
        }
    } 
}
