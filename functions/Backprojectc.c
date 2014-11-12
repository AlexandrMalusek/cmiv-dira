#include <math.h>
#include "mex.h"

/* Input Arguments */
#define P    (prhs[0])
#define THETA  (prhs[1])
#define N_SIZE (prhs[2])
#define INTERP (prhs[3])

/* Output Arguments */
#define  IMG  (plhs[0])

void 
mexFunction(int nlhs, mxArray  *plhs[], int nrhs, const mxArray  *prhs[])
{
  /* INPUT PARAMETERS */
  double *p;                  /* filtered projections (one per col) */
  double *thetaPtr;
  int numAngles;              /* number of projection angles (length of cosines vector, eg) */
  double *Nptr;               /* size of reconstructed image */
  double *interp_ptr;         /* 0 for NN interp, 1 for linear interp */
  
  /* OUTPUT PARAMETERS */
  double *img;                /* output image */

  /* Other variables */
  int N;                      /* integer copy of Nptr (above) */
  int interp_flag;            /* integer copy of interp_ptr (above) */
  int x, y, k;                /* loop indecies */
  double xcoord;              /* stores x-coordinate of pixel */
  double ctr, xleft, ytop;    /* used to tranform from matrix indecies */
                              /* to (x,y) coords (see iradon.m) */
  int len;                    /* length of each projection (spatial dimension) */
  int center;                 /* center index for projections */
  
  /* Temporary variable used for code optimization */
  double cos_theta, sin_theta, t, fraction;
  int a, out_row, angle_row;
  
  /* Check validity of arguments */
  if (nrhs != 4)
  {
    mexErrMsgTxt("Incorrect number of input arguments.");
  }
  if (nlhs != 1)
  {
    mexErrMsgTxt("Incorrect number of output arguments.");
  }
  if (mxIsSparse(P) || mxIsSparse(THETA) || mxIsSparse(N_SIZE))
  {
    mexErrMsgTxt("Sparse inputs not supported.");
  }
  if (!mxIsDouble(P) || !mxIsDouble(THETA) || !mxIsDouble(N_SIZE))
  {
    mexErrMsgTxt("Inputs must be double.");
  }
  
  thetaPtr = mxGetPr(THETA);
  numAngles = mxGetM(THETA) * mxGetN(THETA);
  
  p = mxGetPr(P);
  len = mxGetM(P);
  
  Nptr = mxGetPr(N_SIZE);
  N = (int) *Nptr;
  
  interp_ptr = mxGetPr(INTERP);
  interp_flag = (int) *interp_ptr;
  
  /* Create a matrix for the return argument */
  IMG = mxCreateDoubleMatrix(N, N, mxREAL);
  img = mxGetPr(IMG);  /* Get pointer to the data array */
  
  ctr = floor((N-1) / 2);
  
  xleft = -ctr;
  ytop = ctr;

  center = (int)floor(len/2);  /* center index for projections */ 

  for (k=0;k<numAngles;k++)
  {
    cos_theta = cos(thetaPtr[k]);
    sin_theta = sin(thetaPtr[k]);
    angle_row = k*len;
    
    for (x=0;x<N;x++)
    {
      xcoord = xleft + x;  /* x-coord */
      t = xcoord*cos_theta + ytop*sin_theta;  /* After this, t can simply be decremented by sin_theta each y-iter */
      out_row = x*N;
    
      for (y=0;y<N;y++)
      {      
        a  = ((int) (t + N)) - N;  /* Shifts t to positive values, to avoid using floor */  
        fraction = t - a;
        a +=center;

        img[out_row + y] += fraction*(p[angle_row + a + 1] - p[angle_row + a]) + p[angle_row + a];
      
        t -= sin_theta;  /* decrement t */
      }
    }
  }
}


