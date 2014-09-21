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
MD3(double *Wei3Ptr, double *atte1matPtr, double *atte2matPtr,
    double *att3Ptr, double *dens3Ptr, bool *maskPtr, int isspecial,
    int image_rows, int image_columns, int image_size);

static char rcs_id[] = "$Revision: 1.10 $";

/* Input Arguments */
#define ATTE1MAT  (prhs[0])
#define ATTE2MAT  (prhs[1])
#define ATT3      (prhs[2])
#define DENS3     (prhs[3])
#define MASK      (prhs[4])
#define ISSPECIAL (prhs[5])

/* Output Arguments */
#define WEI3      (plhs[0])

/**
 * Need 6 input arguments: AttE1mat, AttE2mat, Att3, Dens3, mask, isSpecial
 * Produces 1 output: Wei3
 */
void 
mexFunction(int nlhs, mxArray  *plhs[], int nrhs, const mxArray  *prhs[])
{
  /* Input pointers */
  double *atte1matPtr;  /* 511x511 pointer to array of measured attenuation coefficients at energy E1 */
  double *atte2matPtr;  /* 511x511 pointer to array of measured attenuation coefficients at energy E2 */
  double *att3Ptr;      /* 2x2 pointer to array of tabled attenuation coefficients*/
  double *dens3Ptr;     /* 1x2 pointer to mass density of doublet base materials */
  bool *maskPtr;        /* 511x511 pointer to mask defining the tissue to be decomposed*/
  
  
  /* temp pointers for copying values */
  double *from_doublePtr;   
  double *to_doublePtr;   
  int *intPtr;
  bool *from_boolPtr;
  bool *to_boolPtr;
  
  int image_columns;      /* columns of the image to process */
  int image_rows;         /* rows of the image to process */
  int image_size;         /* Total image size, columns * rows */
  
  /* Size of the 2x2 attenuation coefficients */
  int att3_columns;
  int att3_rows;
  int att3_size;
  
  /* Size of the 2x2 density coefficients */
  int dens3_columns;
  int dens3_rows;
  int dens3_size;

  int isspecial;
  
  int k;            /* Loop counter */
  int ndim;         /* Number of dimensions for WEI3 output array */
  int dims[3];      /* Array specifying the dimensions for WEI3 array */
  
  /* Check validity of arguments */
  if (nrhs != 6)
    {
      mexErrMsgTxt("Incorrect number of INPUT arguments.");
    }
  
  if (nlhs != 1)
    {
      mexErrMsgTxt("Incorrect number of OUTPUT arguments.");
    }
  
  if (mxIsSparse(ATTE1MAT) || mxIsSparse(ATTE2MAT) || mxIsSparse(ATT3) || mxIsSparse(DENS3) || mxIsSparse(MASK))
    {
      mexErrMsgTxt("Sparse inputs not supported.");
    }
  
  if (!mxIsDouble(ATTE1MAT) || !mxIsDouble(ATTE2MAT) || !mxIsDouble(ATT3) 
      || !mxIsDouble(DENS3) || !mxIsDouble(ISSPECIAL))
    {
      mexErrMsgTxt("Input must be double.");
    }
  
  if (!mxIsLogical(MASK) )
    {
      mexErrMsgTxt("Mask must be logical.");
    }

  
  /**
   * Array allocation
   */
  
  /* Get the size of the image */
  image_columns = mxGetN(ATTE1MAT);
  image_rows = mxGetM(ATTE1MAT);
  image_size = image_rows * image_columns; 
  
  /* Allocate memory */
  atte1matPtr = (double *) mxCalloc(image_size, sizeof(double));
  atte2matPtr = (double *) mxCalloc(image_size, sizeof(double));
  maskPtr = (bool *) mxCalloc(image_size, sizeof(bool));
  
  att3_columns = mxGetN(ATT3);
  att3_rows = mxGetM(ATT3);
  att3_size = att3_rows * att3_columns;
  att3Ptr = (double *) mxCalloc(att3_size, sizeof(double));
  
  dens3_columns = mxGetN(DENS3);
  dens3_rows = mxGetM(DENS3);
  dens3_size = dens3_rows* dens3_columns; 
  dens3Ptr = (double *) mxCalloc(dens3_size, sizeof(double));
  
  
  /* Copy all data */
  from_doublePtr = mxGetPr(ATTE1MAT);
  to_doublePtr = atte1matPtr;
  for (k = 0; k < image_size; k++)
    *(to_doublePtr++) = *(from_doublePtr++);
  
  from_doublePtr = mxGetPr(ATTE2MAT);
  to_doublePtr = atte2matPtr;
  for (k = 0; k < image_size; k++)
    *(to_doublePtr++) = *(from_doublePtr++);
  
  from_doublePtr = mxGetPr(ATT3);
  to_doublePtr = att3Ptr;
  for (k = 0; k < att3_size; k++)
    *(to_doublePtr++) = *(from_doublePtr++);
  
  from_doublePtr = mxGetPr(DENS3);
  to_doublePtr = dens3Ptr;
  for (k = 0; k < dens3_size; k++)
    *(to_doublePtr++) = *(from_doublePtr++);
  
  from_boolPtr = mxGetLogicals(MASK);
  to_boolPtr = maskPtr;
  for (k = 0; k < image_size; k++)
    *(to_boolPtr++) = *(from_boolPtr++);
    
  /* isSpecial */
  from_doublePtr = mxGetPr(ISSPECIAL);
  isspecial = (int)*from_doublePtr;  
  
  ndim = 3;
  dims[0] = image_rows;
  dims[1] = image_columns;
  dims[2] = 3;
  WEI3 = mxCreateNumericArray(ndim, dims, mxDOUBLE_CLASS, mxREAL);
  
  MD3(mxGetPr(WEI3), atte1matPtr, atte2matPtr, att3Ptr, dens3Ptr,
          maskPtr, isspecial, image_rows, image_columns, image_size);
}

static void 
MD3(double *Wei3Ptr, double *atte1matPtr, double *atte2matPtr,
    double *att3Ptr, double *dens3Ptr, bool *maskPtr, int isspecial,
    int image_rows, int image_columns, int image_size)
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
	  if(*(maskPtr+k*image_rows+l) == 1) 
            {
	      if((isspecial == 0) ||
		 ((*(atte1matPtr+k*image_rows+l) >= *(att3Ptr)) && 
		  (*(atte2matPtr+k*image_rows+l) >= *(att3Ptr + 1))))
                {
		  M[0][0] = ((*(atte1matPtr+k*image_rows+l)) - (*(att3Ptr+0))) / (*(dens3Ptr+0)) -
		    ((*(atte1matPtr+k*image_rows+l)) - (*(att3Ptr+4))) / (*(dens3Ptr+2));
		  
		  M[0][1] = ((*(atte1matPtr+k*image_rows+l)) - (*(att3Ptr+2))) / (*(dens3Ptr+1)) -
		    ((*(atte1matPtr+k*image_rows+l)) - (*(att3Ptr+4))) / (*(dens3Ptr+2));
		  
		  M[1][0] = ((*(atte2matPtr+k*image_rows+l)) - (*(att3Ptr+1))) / (*dens3Ptr) -
		    ((*(atte2matPtr+k*image_rows+l)) - (*(att3Ptr+5))) / (*(dens3Ptr+2));
		  
		  M[1][1] = ((*(atte2matPtr+k*image_rows+l)) - (*(att3Ptr+3))) / (*(dens3Ptr+1)) -
		    ((*(atte2matPtr+k*image_rows+l)) - (*(att3Ptr+5)))/(*(dens3Ptr+2));                    
		  
		  b[0] = -((*(atte1matPtr+k*image_rows+l)) - *(att3Ptr+4)) / (*(dens3Ptr+2));
		  b[1] = -((*(atte2matPtr+k*image_rows+l)) - *(att3Ptr+5)) / (*(dens3Ptr+2));
		  
		  w[0] = 0;
		  w[1] = 0;
		  
		  /* Use Gaussian elimination to solve the linear equation */
		  quota = M[1][0] / M[0][0];
		  M[1][1] = M[1][1] - (M[0][1] * quota);
		  b[1] = b[1] - (b[0] * quota);
		  w[1] = b[1] / (M[1][1]);
		  b[0] = b[0] - (M[0][1] * w[1]);
		  w[0] = b[0] / (M[0][0]);
		  
		  *(Wei3Ptr+k*image_rows+l) = w[0];
		  *(Wei3Ptr+k*image_rows+l+1*image_size) = w[1];        
		  *(Wei3Ptr+k*image_rows+l+2*image_size) = 1.0 - w[0] - w[1]; 
                }
	      else
                {
		  *(Wei3Ptr+k*image_rows+l) = ((*(atte1matPtr+k*image_rows+l) / (*att3Ptr))
					       + (*(atte2matPtr+k*image_rows+l) / (*(att3Ptr+1)))) / 2.0;
                }
            }
        }
    } 
}
