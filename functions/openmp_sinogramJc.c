/*-------------------------------------------------------------------*/
/* This routine computes a Joseph sinogram from a pixelized phantom. */
/* In a Joseph sinogram, the interpolation kernel is designed        */
/* relative to the phantom grid.                                     */
/* Different interpolation filters are available.                    */
/* The routine is based on Matlabs iradon.m.                         */
/*                                                                   */
/* Written by Maria Magnusson Seger 2003-04                          */
/* Updated by Alexander Örtenberg   2014-11                          */
/*-------------------------------------------------------------------*/
#include <math.h>
#include <omp.h>
#include "mex.h"

static void sinogramJ(double *pPtr, double *iPtr, double *thetaPtr, double *rinPtr,
          int M, int N, int xOrigin, int yOrigin, int numAngles, int rFirst, 
          int rSize, int interpolation);

static char rcs_id[] = "$Revision: 1.10 $";

#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define MIN(x,y) ((x) < (y) ? (x) : (y))

#define PI 3.14159265358979

/* Input Arguments */
#define I      (prhs[0])
#define THETA  (prhs[1])
#define R_IN   (prhs[2])
#define INTERP (prhs[3])

/* Output Arguments */
#define  P      (plhs[0])
#define R      (plhs[1])

void 
mexFunction(int nlhs, mxArray  *plhs[], int nrhs, const mxArray  *prhs[])
{
  int numAngles;    /* number of theta values */
  int numProjval;    /* number of projection values */
  double *thetaPtr;    /* pointer to theta values in radians */
  double *rinPtr;    /* pointer to projection coordinate array */
  double *pr1, *pr2;  /* help pointers used in loop */
  double deg2rad;    /* conversion factor */
  int k;                /* loop counter */
  int M, N;             /* input image size */
  int xOrigin, yOrigin;  /* center of image */
  int rFirst, rLast;  /* r-values for first and last row of output */
  int rSize;      /* number of rows in output */
  int interpolation;  /* interpolation type */

  /* Check validity of arguments */
  if (nrhs < 4)
  {
      mexErrMsgTxt("Too few input arguments");
  }
  if (nrhs > 4)
  {
      mexErrMsgTxt("Too many input arguments");
  }
  if (nlhs > 2)
  {
      mexErrMsgTxt("Too many output arguments to SINOGRAMJ");
  }
  if (mxIsSparse(I) || mxIsSparse(THETA) || mxIsSparse(R_IN))
  {
      mexErrMsgTxt("Sparse inputs not supported");
  }
  if (!mxIsDouble(I) || !mxIsDouble(THETA) || !mxIsDouble(R_IN) || !mxIsDouble(INTERP))
  {
      mexErrMsgTxt("Inputs must be double");
  }

  /* Get THETA degree values and convert to radians */
  deg2rad = PI / 180.0;
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

  /* Get INTERP values */
  pr1 = mxGetPr(INTERP);
  interpolation = *pr1;

  /* Get input image size */
  M = mxGetM(I);
  N = mxGetN(I);

  /* Where is the coordinate system's origin? */
  xOrigin = MAX(0, (N-1)/2);
  yOrigin = MAX(0, (M-1)/2);

  /* Second out parameter? */
  if (nlhs == 2)
  {
    R = mxCreateDoubleMatrix(rSize, 1, mxREAL);
    pr1 = mxGetPr(R);
    for (k = rFirst; k <= rLast; k++)
      *(pr1++) = (double) k;
  }
  
  /* Invoke main computation routines */
  if (mxIsComplex(I))
  {
    P = mxCreateDoubleMatrix(rSize, numAngles, mxCOMPLEX);
    sinogramJ(mxGetPr(P), mxGetPr(I), thetaPtr, rinPtr, M, N, xOrigin, yOrigin, 
       numAngles, rFirst, rSize, interpolation); 
    sinogramJ(mxGetPi(P), mxGetPi(I), thetaPtr, rinPtr, M, N, xOrigin, yOrigin, 
       numAngles, rFirst, rSize, interpolation);
  }
  else
  {
    P = mxCreateDoubleMatrix(rSize, numAngles, mxREAL);
    sinogramJ(mxGetPr(P), mxGetPr(I), thetaPtr, rinPtr, M, N, xOrigin, yOrigin, 
       numAngles, rFirst, rSize, interpolation);
  }
}

static void 
sinogramJ(double *pPtr, double *iPtr, double *thetaPtr, double *rinPtr, int M, int N, 
    int xOrigin, int yOrigin, int numAngles, int rFirst, int rSize, int interpolation)
{
    
  int x,y,k,i;    /* Loop variables */
  int radius;     /* Radius of circle from which to use values */
  double r;     /* Polar coordinate */
  int r_index;    /* Polar coordinate as integer to index matrix */
  double fraction;  /* Fraction of the r coordinate */
  
  double pixelvalue, slopedpixelvalue;  /* Pixelvalue from input image and scaled version */
  double leftpixel, rightpixel;     /* Distribution for left and right pixel */
  double distance, leftdistance, rightdistance; /* Distance to left and right pixel */
  
  int *xdistance, *ydistance;/*Distance in carthesian coordinates to image center */
  int *pixelindices;         /* Store indices of pixels to calculate */
  int xdist, ydist; /* temporary variables */
  int pixelindex; /* Current index to store pixel data on */
  double pixelradius; /* Radius of the pixel from center of the image */

  /* Precalculate the values for all angles */
  double angle;
  double *cosine, *sine, *slope;
  
  cosine  = (double *) malloc(numAngles * sizeof(double));
  sine    = (double *) malloc(numAngles * sizeof(double));
  slope   = (double *) malloc(numAngles * sizeof(double));
  
  for(k=0;k< numAngles;++k)
  {
    angle    = -thetaPtr[k];
    cosine[k] = cos(angle);
    sine[k]   = sin(angle);
    /* Calculate the slope depending on which angle value is larger */
    slope[k]   = 1/MAX(fabs(cosine[k]), fabs(sine[k]));
  }
  
  /* Only values in a circle will be used, the edges do not add anything */
  radius = ceil(rSize/2);  
  
  xdistance    = (int *)malloc (sizeof(int) * M * N);
  ydistance    = (int *)malloc (sizeof(int) * M * N);
  pixelindices = (int *)malloc (sizeof(int) * M * N);
  pixelindex   = 0;
  
  /** Checks for every pixel if it is within the radius of the unit circle
   *  and if it is a non-zero value. Only store its index and values if it
   *  passes both checks, or it does not contribute
   */

  for(y=0;y<M;++y)
  {    
    for(x=0;x<N;++x)
    {
      ydist = x - xOrigin;
      xdist = y - yOrigin;
      pixelradius = sqrt(xdist * xdist + ydist * ydist);
      
      if((iPtr[y*M + x] != 0) && (pixelradius <= radius))
      {
        ydistance[pixelindex] = ydist;
        xdistance[pixelindex] = xdist;
        pixelindices[pixelindex] = y*N+x;
        ++pixelindex;
      }
    }
  }
  
  #pragma omp parallel for private(i,pixelvalue, r, r_index, fraction, distance,\
                                   leftdistance, rightdistance, slopedpixelvalue,\
                                   leftpixel, rightpixel)
  /* Calculate for every angle given as input*/
  for(k=0;k<numAngles;++k)
  {
    /* Calculate for all pixels that will contribute */
    for(i=0;i<pixelindex;++i)
    {
      /* Inside the circle, get the pixel value */
      pixelvalue = iPtr[pixelindices[i]];
          
      /* Find the index for the radial coordinates */
      r = xdistance[i]*cosine[k] + ydistance[i]*sine[k];          
      r += xOrigin;   /* add xOrigin to shift center of image, avoiding negative values */
      r_index = (int) r;  
      fraction = r - r_index;

      /* Get the pixel value and distribute between two pixels
       * Calculates the distance once as it is used multiple times
       * The slope is dependent on the angle, decreasing the
       * triangle size */
      distance = fraction*slope[k];
      /* No contribution if the distance is less than 0 */
      /* Equal to 
       * (1 - fraction*slope[k]) and
       * (1 - (1 - fraction) * slope[k])*/
      leftdistance  = MAX(0, (1 - distance));
      rightdistance = MAX(0, (1 + distance - slope[k]));
  
      slopedpixelvalue = pixelvalue * slope[k];
      leftpixel  = leftdistance  * slopedpixelvalue;
      rightpixel = rightdistance * slopedpixelvalue;

      pPtr[k*M + r_index] += leftpixel;
      pPtr[k*M + r_index + 1] += rightpixel;
    }
  }

  free(xdistance);
  free(ydistance);
  free(pixelindices);
  free(cosine);
  free(sine);
  free(slope);
}

  