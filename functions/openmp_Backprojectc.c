#include <math.h>
#include <omp.h>
#include "mex.h"

/* Input Arguments */
#define P      (prhs[0])
#define THETA  (prhs[1])
#define N_SIZE (prhs[2])
#define INTERP (prhs[3])

/* Output Arguments */
#define	IMG    (plhs[0])

void 
mexFunction(int nlhs, mxArray  *plhs[], int nrhs, const mxArray  *prhs[])
{
	/* INPUT PARAMETERS */
  double *p;                 /* filtered projections (one per col) */
  double *thetaPtr;
  int numAngles;             /* number of projection angles (length of cosines vector, eg) */
  double angle;
	double *cosines;           /* pre-computed cos of proj angles */
	double *sines;             /* pre-computed sin of proj angles */
	double *Nptr;              /* size of reconstructed image */
	double *interp_ptr;        /* 0 for NN interp, 1 for linear interp */
	
	/* OUTPUT PARAMETERS */
	double *img;               /* output image */

	/* Other variables */
	int N;                     /* integer copy of Nptr (above) */
	int interp_flag;           /* integer copy of interp_ptr (above) */
	int x, y, i, k;            /* loop indecies */
	double xcoord;             /* stores x-coordinate of pixel */
	double ctr, xleft, ytop;	 /* used to tranform from matrix indecies */
                             /* to (x,y) coords (see iradon.m) */
	int len;                   /* length of each projection (spatial dimension) */
	int ctr_idx;               /* centre index for projections */
	
	/* Temporary variable used for code optimization */
	double cos_theta, sin_theta, t;
	int a;
	double *proj;              /* points at the start of a projection (a column) */
	double *private_img;
  double *temp_img;
  int thread_count;
  
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
  
  cosines = (double *) malloc(numAngles * sizeof(double));
  sines   = (double *) malloc(numAngles * sizeof(double));
  
  /* Precalculate all cosine and sine values */
  for(k=0;k< numAngles;++k)
  {
    angle = thetaPtr[k];
    cosines[k] = cos(angle);    /* Calculate cosine value */
    sines[k] = sin(angle);      /* Calculate sine value */
  }
	
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
	
	ctr_idx = (int)floor(len/2);  /* centre index for projections */ 
  
  thread_count = omp_get_max_threads();
  
  temp_img = (double *) malloc(sizeof(double)* N * N * thread_count);
  for(i=0;i<N*N*thread_count;++i)
    temp_img[i] = 0;
  
  for(i=0;i<N*N;++i)
    img[i] = 0;
  
  #pragma omp parallel private(private_img)
  {     
    private_img = temp_img + N*N*omp_get_thread_num();
    double *img_ptr;
	
    #pragma omp for private(k, cos_theta, sin_theta, proj, img_ptr, x, xcoord, t, y, a)
    for (k=0;k<numAngles;k++)
    {
      cos_theta = cosines[k];
      sin_theta = sines[k];
      proj = (p + k*len);  /* point at proper column */
      img_ptr = private_img;

      for (x=0;x<N;x++)
      {
        xcoord = xleft + x;  /* x-coord */
        t = xcoord*cos_theta + ytop*sin_theta;  /* After this, t can simply be decremented by sin_theta each y-iter */

        for (y=0;y<N;y++)
        {			
          a = ((int) (t + N)) - N;  /* Shifts t to positive values, to avoid using floor */	

          *img_ptr += (t-a)*( proj[a + ctr_idx + 1] - proj[a+ctr_idx] ) + proj[a+ctr_idx];
          img_ptr++;

          t -= sin_theta;
        }
      }
    }
    
    #pragma omp critical
    {
      for(i=0;i<N*N;++i)
      {
        img[i] += private_img[i];
      }
    }
  }
  
  free(temp_img);
  free(cosines);
  free(sines);
}