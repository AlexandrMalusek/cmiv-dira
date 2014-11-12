#include <math.h>
#include "mex.h"

static void WBHC(double *rebsimPtr, double *polycrPtr, double *distPtr, int M, int N, int polysize);

static char rcs_id[] = "$Revision: 1.10 $";

/* Input Arguments */
#define N1     (prhs[0])
#define M1     (prhs[1])
#define REBSIM (prhs[2])
#define POLYCR (prhs[3])

/* Output Arguments */
#define	DIST   (plhs[0])

void 
mexFunction(int nlhs, mxArray  *plhs[], int nrhs, const mxArray  *prhs[])
{
  double *rebsimPtr;		/* pointer to theta values in radians */
  double *polycrPtr;		/* pointer to projection coordinate array */
  double *pr1, *pr2;	/* help pointers used in loop */
  int k, size;                /* loop counter */
  int M, N, polysize;             /* input image size */
  


  /* Check validity of arguments */
  if (nrhs != 4)
  {
      mexErrMsgTxt("Incorrect number of input arguments");
  }
  if (nlhs != 1)
  {
      mexErrMsgTxt("Incorrect number of output arguments");
  }
  if (mxIsSparse(M1) || mxIsSparse(N1) || mxIsSparse(REBSIM) || mxIsSparse(POLYCR))
  {
      mexErrMsgTxt("Sparse inputs not supported");
  }
  if (!mxIsDouble(M1) || !mxIsDouble(N1) || !mxIsDouble(REBSIM) || !mxIsDouble(POLYCR))
  {
      mexErrMsgTxt("Inputs must be double");
  }
  
  pr1 = mxGetPr(M1);
  M = (int) *pr1;
  
  pr1 = mxGetPr(N1);
  N = (int) *pr1;
  
  polysize = mxGetM(POLYCR) * mxGetN(POLYCR);
  polycrPtr = (double *) mxCalloc(polysize, sizeof(double));
  pr1 = mxGetPr(POLYCR);
  pr2 = polycrPtr;
  for(k=0;k<polysize;++k)
    *(pr2++) = *(pr1++);
    
  
  size = mxGetM(REBSIM) * mxGetN(REBSIM);
  rebsimPtr = (double *) mxCalloc(size, sizeof(double));
  pr1 = mxGetPr(REBSIM);
  pr2 = rebsimPtr;
  for(k=0;k<size;++k)
    *(pr2++) = *(pr1++);
    
  DIST = mxCreateDoubleMatrix(N, M, mxREAL);
  
  WBHC(rebsimPtr, polycrPtr, mxGetPr(DIST), M, N, polysize);
}

static void WBHC(double *rebsimPtr, double *polycrPtr, double *distPtr, int M, int N, int polysize)
{
    int i, j, k;
    double value;

    for(i=0;i<M;++i)
    {
        for(j=0;j<N;++j)
        {
            value = rebsimPtr[i*N + j];
            
            for(k=0;k<polysize-2;++k)
            {
                
                if(k == 0 && value < polycrPtr[k])
                {
                    distPtr[i*N + j] = value / polycrPtr[0];
                    break;
                }
                else if(( polycrPtr[k] <= value) && (value < polycrPtr[k+1]))
                {
                    distPtr[i*N + j] = (k + 1) + (value - polycrPtr[k]) / (polycrPtr[k+1] - polycrPtr[k]); 
                    break;
                }
            }
        }
    }
}

	
