/*
 * Backprojection is performed by placing the projections values
 * along a line integral (radiological path). This is performed
 * for all projections, for all angles. This implementation works
 * by iterating over all projections and placing the values on a single
 * output row and then moving on to the next output row. The benefit of
 * this implementation is that for parallel implementations no
 * synchronization is required and no temporary data needs to be allocated.
 * This work is based on the implementation made by Jeff Orchard, 
 * http://www.mathworks.com/matlabcentral/fileexchange/12852-iradon-speedy
 *
 * Created by Alexander Örtenberg 2015-04
 */
#include <math.h>
#include "mex.h"

/* Input Arguments */
#define P      (prhs[0])
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
  double *thetaPtr;           /* pointer to array of projection angles */
  int numAngles;              /* number of projection angles (length of cosines vector, eg) */
  double *Nptr;               /* size of reconstructed image */
  double *interp_ptr;         /* 0 for NN interp, 1 for linear interp */
  
  double *img;                /* output image */

  int N;                      /* integer copy of Nptr (above) */
  int interp_flag;            /* integer copy of interp_ptr (above) */
  int x, y, k;                /* loop indecies */
  double xcoord;              /* stores x-coordinate of pixel */
  double ctr, xleft, ytop;    /* used to tranform from matrix indecies */
                              /* to (x,y) coords (see iradon.m) */
  int projection_length;      /* length of each projection (spatial dimension) */
  int center;                 /* center index for projections */
  
  /* Temporary variable used for code optimization */
  double cos_theta, sin_theta, t, fraction;
  int a, out_row, input_row;
  
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
  
  /* Read the data and values needed */
  thetaPtr = mxGetPr(THETA);
  numAngles = mxGetM(THETA) * mxGetN(THETA);
  
  p = mxGetPr(P);
  projection_length= mxGetM(P);
  
  Nptr = mxGetPr(N_SIZE);
  N = (int) *Nptr;
  
  interp_ptr = mxGetPr(INTERP);
  interp_flag = (int) *interp_ptr;
  
  /* Create a matrix for the return argument */
  IMG = mxCreateDoubleMatrix(N, N, mxREAL);
  img = mxGetPr(IMG);
  
  ctr = floor((N-1) / 2);
  
  xleft = -ctr;
  ytop = ctr;
  /* center index for projections */ 
  center = (int)floor(projection_length/2);
  
  /* For each row in the output matrix */
  for(x=0;x<N;x++)
  {
    xcoord = xleft + x;
    out_row = x*N;
    
    /* For every projection angle in the input matrix (each input matrix row) */
    for(k=0;k<numAngles;k++)
    {
      cos_theta = cos(thetaPtr[k]);
      sin_theta = sin(thetaPtr[k]);
      
      /* Set the counter to current input row*/
      input_row = k*projection_length;
      
      t = xcoord*cos_theta + ytop*sin_theta; 
      
      for (y=0;y<N;y++)
      {      
        /* Calculate what values to place from the input matrix */
        a  = ((int) (t + N)) - N;  /* Shifts t to positive values, to avoid using floor */  
        fraction = t - a;
        a +=center;

        img[out_row + y] += fraction*(p[input_row + a + 1] - p[input_row + a]) + p[input_row + a];
      
        /* Step forward to next position to read data from */
        t -= sin_theta;
      }
    }
  }
}




