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
#include <omp.h>
#include "mex.h"

static void sinogramJ(double *pPtr, double *iPtr, double *thetaPtr, double *rinPtr,
		      int M, int N, int xOrigin, int yOrigin, int numAngles, int rFirst, 
		      int rSize, int filter);

static char rcs_id[] = "$Revision: 1.10 $";

#define MAXX(x,y) ((x) > (y) ? (x) : (y))

#define Pi 3.14159265358979
#define SMALL 0.000000001

/* Input Arguments */
#define I      (prhs[0])
#define THETA  (prhs[1])
#define R_IN   (prhs[2])
#define FILTER (prhs[3])

/* Output Arguments */
#define	P      (plhs[0])
#define R      (plhs[1])

void 
mexFunction(int nlhs, mxArray  *plhs[], int nrhs, const mxArray  *prhs[])
{
  int numAngles;		/* number of theta values */
  int numProjval;		/* number of projection values */
  double *thetaPtr;		/* pointer to theta values in radians */
  double *rinPtr;		/* pointer to projection coordinate array */
  double *pr1, *pr2;		/* help pointers used in loop */
  double deg2rad;		/* conversion factor */
  int k;			/* loop counter */
  int M, N;			/* input image size */
  int xOrigin, yOrigin;		/* center of image */
  int rFirst, rLast;		/* r-values for first and last row of output */
  int rSize;			/* number of rows in output */
  int filter;			/* filter type */
  
  /* Check validity of arguments */
  if (nrhs < 4) {
    mexErrMsgTxt("Too few input arguments");
  }
  if (nrhs > 4) {
    mexErrMsgTxt("Too many input arguments");
  }
  if (nlhs > 2) {
    mexErrMsgTxt("Too many output arguments to SINOGRAMJ");
  }
  if (mxIsSparse(I) || mxIsSparse(THETA) || mxIsSparse(R_IN)) {
    mexErrMsgTxt("Sparse inputs not supported");
  }
  if (!mxIsDouble(I) || !mxIsDouble(THETA) || !mxIsDouble(R_IN) || !mxIsDouble(FILTER)) {
    mexErrMsgTxt("Inputs must be double");
  }
  
  /* Get THETA values */
  deg2rad = 3.14159265358979 / 180.0;
  numAngles = mxGetM(THETA) * mxGetN(THETA);
  thetaPtr = (double *) mxCalloc(numAngles, sizeof(double));
  pr1 = mxGetPr(THETA);
  pr2 = thetaPtr;
  for (k = 0; k < numAngles; k++)
    *(pr2++) = *(pr1++) * deg2rad;
  
  
  /* Get R_IN values */
  numProjval = mxGetM(R_IN) * mxGetN(R_IN);
  rinPtr = (double *) mxCalloc(numProjval, sizeof(double));
  pr1 = mxGetPr(R_IN);
  pr2 = rinPtr;
  for (k = 0; k < numProjval; k++)
    *(pr2++) = *(pr1++);
  rSize  = numProjval;
  rFirst = (1-rSize)/2;
  rLast  = -rFirst;
  
  /* Get FILTER values */
  pr1 = mxGetPr(FILTER);
  filter = *pr1;

  /* Get input image size */
  M = mxGetM(I);
  N = mxGetN(I);

  /* Where is the coordinate system's origin? */
  xOrigin = MAXX(0, (N-1)/2);
  yOrigin = MAXX(0, (M-1)/2);

  /* Second out parameter? */
  if (nlhs == 2) {
    R = mxCreateDoubleMatrix(rSize, 1, mxREAL);
    pr1 = mxGetPr(R);
    for (k = rFirst; k <= rLast; k++)
      *(pr1++) = (double) k;
  }
  
  /* Invoke main computation routines */
  if (mxIsComplex(I)) {
    P = mxCreateDoubleMatrix(rSize, numAngles, mxCOMPLEX);
    sinogramJ(mxGetPr(P), mxGetPr(I), thetaPtr, rinPtr, M, N, xOrigin, yOrigin, 
	      numAngles, rFirst, rSize, filter); 
    sinogramJ(mxGetPi(P), mxGetPi(I), thetaPtr, rinPtr, M, N, xOrigin, yOrigin, 
	      numAngles, rFirst, rSize, filter);
  } else {
    P = mxCreateDoubleMatrix(rSize, numAngles, mxREAL);
    sinogramJ(mxGetPr(P), mxGetPr(I), thetaPtr, rinPtr, M, N, xOrigin, yOrigin, 
	      numAngles, rFirst, rSize, filter);
  }
}

static void 
sinogramJ(double *pPtr, double *iPtr, double *thetaPtr, double *rinPtr, int M, int N, 
	  int xOrigin, int yOrigin, int numAngles, int rFirst, int rSize, int filter)
{
  int    k, n, m, p;        /* loop counters */
  double k_s, l_s, m_s;	    /* Siddon variables */
  double angle;		    /* radian angle value */
  double cosine, sine;      /* cosine and sine of current angle */
  double *pr;               /* points inside output array */
  double *pixelPtr;         /* points inside input array */
  double r;		    /* radial coordinate */
  double pixel1, pixel2;    /* current pixel values */
  double pixel3, pixel0;    /* current pixel values */
  double mIdx, nIdx;	    /* m and n value offset from initial array element */
  double loc1, loc2;        /* pixel locations to the left and right */
  double loc0, loc3;        /* more locations */
  double value;		    /* help value */
  int    mLow, nLow;	    /* (int) mIdx, (int) nIdx */

  /*printf("filter=%d\n", filter);
    printf("rSize=%d\n", rSize);*/
  
  /*-------------------------------------*/
  /* Projection generation MAIN LOOP.    */
  /* Loop through all projection angles. */
  /*-------------------------------------*/
  
#pragma omp parallel for private(angle, cosine, sine, pr, p,  r, m, n, pixelPtr, mIdx, mLow, loc1, loc2, pixel1, pixel2, nIdx, nLow)
  for (k = 0; k < numAngles; k++)
    {
      angle  = -thetaPtr[k];
      cosine = cos(angle); 
      sine   = sin(angle);
      pr     = pPtr + k*rSize;	/* pointer to the top of the output column */
      /* linear interpolation */
      if (fabs(cosine)<fabs(sine))
	{
	  for (p = 0; p < rSize; p++)
	    {
	      r = rinPtr[p];
	      for (n = 0; n < N; n++)
		{
		  pixelPtr = iPtr + n*M;
		  mIdx = (r - (n-xOrigin) * cosine) / sine + yOrigin;  
		  if ((mIdx>0) && (mIdx<(M-1)))
		    {
		      mLow = (int) mIdx;
		      loc1 = mIdx - mLow;
		      loc2 = 1 - loc1;
		      pixel1 = *(pixelPtr+mLow);
		      pixel2 = *(pixelPtr+mLow+1);
		      pr[p] += (pixel1 * (1-loc1) + pixel2 * (1-loc2)) / fabs(sine);
		    }
		}
	    }
	} 
      else
	{
	  for (p = 0; p < rSize; p++)
	    {
	      r = rinPtr[p];
	      for (m = 0; m < M; m++)
		{
		  pixelPtr = iPtr + m;
		  nIdx = (r - (m-yOrigin) * sine) / cosine + xOrigin;  
		  if ((nIdx>0) && (nIdx<(N-1)))
		    {
		      nLow = (int) nIdx;
		      loc1 = nIdx - nLow;
		      loc2 = 1 - loc1;
		      pixel1 = *(pixelPtr+M*nLow);
		      pixel2 = *(pixelPtr+M*(nLow+1));
		      pr[p] += (pixel1 * (1-loc1) + pixel2 * (1-loc2)) / fabs(cosine);
		    }
		}
	    }
	}    
    }
}
