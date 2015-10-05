/*-------------------------------------------------------------------*/
/* This routine computes a Joseph sinogram from a pixelized phantom. */
/* In a Joseph sinogram, the interpolation kernel is designed        */
/* relative to the phantom grid.                                     */
/*                                                                   */
/* Written by Maria Magnusson Seger 2003-04                          */
/* Updated by Alexander Örtenberg   2014-11                          */
/*-------------------------------------------------------------------*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <CL/cl.h>
#include "mex.h"

static void sinogramJ(double *pPtr, double *iPtr, double *thetaPtr, double *rinPtr,
          int M, int N, int xOrigin, int yOrigin, int numAngles, int rFirst, 
          int rSize, int interpolation);

static char rcs_id[] = "$Revision: 1.10 $";

#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define PI 3.14159265358979

/* Input Arguments */
#define I      (prhs[0])
#define THETA  (prhs[1])
#define R_IN   (prhs[2])
#define INTERP (prhs[3])

/* Output Arguments */
#define  P      (plhs[0])
#define  R      (plhs[1])

void 
mexFunction(int nlhs, mxArray  *plhs[], int nrhs, const mxArray  *prhs[])
{
  int numAngles;        /* number of theta values */
  int numProjval;       /* number of projection values */
  double *thetaPtr;     /* pointer to theta values in radians */
  double *rinPtr;       /* pointer to projection coordinate array */
  double *pr1, *pr2;    /* help pointers used in loop */
  double deg2rad;       /* conversion factor */
  int k;                /* loop counter */
  int M, N;             /* input image size */
  int xOrigin, yOrigin; /* center of image */
  int rFirst, rLast;    /* r-values for first and last row of output */
  int rSize;            /* number of rows in output */
  int interpolation;    /* interpolation type */

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
  int x,y,k;                                     /* Loop variables */
  int radius;                                    /* Radius of circle from which to use values */
  
  /* Pixel information */
  int *xdistance, *ydistance;                    /*Distance in carthesian coordinates to image center */
  int *pixelindices;                             /* Store indices of pixels to calculate */
  int xdist, ydist;                              /* temporary variables */
  int pixelindex;                                /* Current index to store pixel data on */
  double pixelradius;                            /* Radius of the pixel from center of the image */
    
  /* Precalculate the values for all angles */
  double angle;
  double *cosine, *sine, *slope;
  
  cosine  = (double *) malloc(numAngles * sizeof(double));
  sine    = (double *) malloc(numAngles * sizeof(double));
  slope   = (double *) malloc(numAngles * sizeof(double));
  
  for(k=0;k< numAngles;++k)
  {
    angle     = -thetaPtr[k];
    cosine[k] = cos(angle);
    sine[k]   = sin(angle);
    /* Calculate the slope depending on which angle value is larger */
    slope[k]   = 1/MAX(fabs(cosine[k]), fabs(sine[k]));
  }
  
  /* Only values in a circle will be used, the edges do not add anything */
  radius = rSize/2;  
  
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
      pixelradius = sqrt((double) (xdist * xdist + ydist * ydist));
      
      if((iPtr[y*M + x] != 0) && (pixelradius <= radius))
      {
        ydistance[pixelindex] = ydist;
        xdistance[pixelindex] = xdist;
        pixelindices[pixelindex] = y*N+x;
        ++pixelindex;
      }
    }
  } 
  
  /* OpenCL variables */
  int error;
  cl_platform_id platform;
  unsigned int no_platforms;
  cl_device_id device;
  unsigned int no_devices;
  size_t no_workgroups;
  cl_context context;
  cl_command_queue commandqueue;
  
  static cl_program program_sinogramJ;
  static cl_kernel kernel_sinogramJ;
  
  cl_mem IMG_input; 
  cl_mem IDX_input;
  cl_mem XDIS_input;
  cl_mem YDIS_input;
  cl_mem COS_input;
  cl_mem SIN_input;
  cl_mem SLOPE_input;
  cl_mem PROJ_output;

  size_t localworksize;
  size_t globalworksize;
  
  cl_event event;
  
  /* Kernel to execute */
  const char *source =   "\n" \
"#define MAX(x,y) ((x) > (y) ? (x) : (y))  \n" \
"#define MIN(x,y) ((x) < (y) ? (x) : (y))  \n" \
"  \n" \
"__kernel void sinogramJ(const int no_pixels, __global double *image,  \n" \
"                        __global int *indices, __global int *xdistance,  \n" \
"                        __global int *ydistance, __global double *cosine,  \n" \
"                        __global double *sine, const int xOrigin,  \n" \
"                        __global double *slope, const int M, __global double *projections)  \n" \
"{   \n" \
"  int i;  \n" \
"  double pixelvalue;  \n" \
"  double r, fraction;  \n" \
"  int r_index;  \n" \
"  double distance;  \n" \
"  double leftdistance, rightdistance;  \n" \
"  double slopedpixelvalue, leftpixel, rightpixel;  \n" \
"  \n" \
"  double max_r = 0;  \n" \
"  \n" \
"  int k = get_global_id(0);  \n" \
"            \n" \
"  for(i=0;i<M;++i)  \n" \
"  {  \n" \
"    projections[k*M + i] = 0;  \n" \
"  }  \n" \
"  \n" \
"  for(i=0;i<no_pixels;++i)  \n" \
"  {  \n" \
"    pixelvalue = image[indices[i]];  \n" \
"  \n" \
"    r = xdistance[i]*cosine[k] + ydistance[i]*sine[k];        \n" \
"  \n" \
"    r += xOrigin;     \n" \
"    r_index = (int) r;    \n" \
"    fraction = r - r_index;  \n" \
"  \n" \
"    distance = fraction*slope[k];  \n" \
"  \n" \
"    leftdistance  = MAX(0, (1 - distance));  \n" \
"    rightdistance = MAX(0, (1 + distance - slope[k]));  \n" \
"  \n" \
"    slopedpixelvalue = pixelvalue * slope[k];  \n" \
"    leftpixel  = leftdistance  * slopedpixelvalue;  \n" \
"    rightpixel = rightdistance * slopedpixelvalue;  \n" \
"  \n" \
"    projections[k*M + r_index + 0] += leftpixel;  \n" \
"    projections[k*M + r_index + 1] += rightpixel;  \n" \
"  \n" \
"  }  \n" \
"}  \n" \
"\n";

  /* Initialization */
  error = clGetPlatformIDs(1, &platform, &no_platforms);
  error = clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 1, &device, &no_devices); 
  error = clGetDeviceInfo(device, CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(size_t), &no_workgroups, NULL);

  context = clCreateContext(0, 1, &device, NULL, NULL, &error);
  commandqueue = clCreateCommandQueue(context, device, 0, &error);
  
  program_sinogramJ = clCreateProgramWithSource(context, 1, (const char **)&source, 
                                                    NULL, &error);
    
  error = clBuildProgram(program_sinogramJ, 0, NULL, NULL, NULL, NULL);
  
  kernel_sinogramJ = clCreateKernel(program_sinogramJ, "sinogramJ", &error);
  
  /* Create buffers for data */
  IMG_input   = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,  sizeof(double) * M * N, iPtr, NULL);
  IDX_input   = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,  sizeof(int) * M * N, pixelindices, NULL);
  XDIS_input  = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,  sizeof(int) * M * N, xdistance, NULL);
  YDIS_input  = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,  sizeof(int) * M * N, ydistance, NULL);
  COS_input   = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,  sizeof(double) * numAngles, cosine, NULL);
  SIN_input   = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,  sizeof(double) * numAngles, sine, NULL);
  SLOPE_input = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,  sizeof(double) * numAngles, slope, NULL);
  PROJ_output = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(double) * rSize * numAngles, NULL, NULL);
  
  /* Set the kernel arguments */
  error  = clSetKernelArg(kernel_sinogramJ, 0, sizeof(int), &pixelindex);
  error |= clSetKernelArg(kernel_sinogramJ, 1, sizeof(cl_mem), &IMG_input);
  error |= clSetKernelArg(kernel_sinogramJ, 2, sizeof(cl_mem), &IDX_input);
  error |= clSetKernelArg(kernel_sinogramJ, 3, sizeof(cl_mem), &XDIS_input);
  error |= clSetKernelArg(kernel_sinogramJ, 4, sizeof(cl_mem), &YDIS_input);
  error |= clSetKernelArg(kernel_sinogramJ, 5, sizeof(cl_mem), &COS_input);
  error |= clSetKernelArg(kernel_sinogramJ, 6, sizeof(cl_mem), &SIN_input);
  error |= clSetKernelArg(kernel_sinogramJ, 7, sizeof(int), &xOrigin);
  error |= clSetKernelArg(kernel_sinogramJ, 8, sizeof(cl_mem), &SLOPE_input);
  error |= clSetKernelArg(kernel_sinogramJ, 9, sizeof(int), &rSize);
  error |= clSetKernelArg(kernel_sinogramJ, 10, sizeof(cl_mem), &PROJ_output);
  mexPrintf("arguments error: %i \n", error);
  
  localworksize = 32;
  globalworksize = numAngles;

  /* Enqueue and run */
  error = clEnqueueNDRangeKernel(commandqueue, kernel_sinogramJ, 1, NULL, &globalworksize, &localworksize, 0, NULL, &event);
  clWaitForEvents(1, &event); // Synch
  
    /* Read result */
  error = clEnqueueReadBuffer(commandqueue, PROJ_output, CL_TRUE, 0, sizeof(double) * rSize * numAngles, pPtr, 0, NULL, NULL);
  clWaitForEvents(1, &event); // Synch
  
  clReleaseMemObject(IMG_input);
  clReleaseMemObject(IDX_input);
  clReleaseMemObject(XDIS_input);
  clReleaseMemObject(YDIS_input);
  clReleaseMemObject(COS_input);
  clReleaseMemObject(SIN_input);
  clReleaseMemObject(SLOPE_input);
  clReleaseMemObject(PROJ_output);
  
  clReleaseKernel(kernel_sinogramJ);
  clReleaseProgram(program_sinogramJ);
  clReleaseCommandQueue(commandqueue);
  clReleaseContext(context);
  
  free(xdistance);
  free(ydistance);
  free(pixelindices);
  free(cosine);
  free(sine);
  free(slope);
}
