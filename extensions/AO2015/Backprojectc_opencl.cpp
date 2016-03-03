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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <CL/cl.h>
#include "mex.h"

/* Input Arguments */
#define P      (prhs[0])
#define THETA  (prhs[1])
#define N_SIZE (prhs[2])
#define INTERP (prhs[3]) /* Not currently used */

/* Output Arguments */
#define  IMG    (plhs[0])

void 
mexFunction(int nlhs, mxArray  *plhs[], int nrhs, const mxArray  *prhs[])
{
  /* INPUT PARAMETERS */
  double *projections;        /* filtered projections (one per col) */
  double *thetaPtr;           /* pointer to array of projection angles */
  int numAngles;              /* number of projection angles (length of cosines vector, eg) */
  double *Nptr;               /* size of reconstructed image */
  double *interp_ptr;         /* 0 for NN interp, 1 for linear interp */
  
  double *imgPtr;             /* output image */

  int N;                      /* integer copy of Nptr (above) */
  int interp_flag;            /* integer copy of interp_ptr (above) */
  int x, y, k;                /* loop indecies */
  double xcoord;              /* stores x-coordinate of pixel */
  double center, xleft, ytop; /* used to tranform from matrix indecies */
                              /* to (x,y) coords (see iradon.m) */
  int projection_length;      /* length of each projection (spatial dimension) */
  
  /* Temporary variable used for code optimization */
  double cos_theta, sin_theta, t, fraction;
  
  /* OpenCL variables */
  int error;
  cl_platform_id platform;
  unsigned int no_platforms;
  cl_device_id device;
  unsigned int no_devices;
  size_t no_workgroups;
  cl_context context;
  cl_command_queue commandqueue;
  
  static cl_program program_Backproject;
  static cl_kernel kernel_Backproject;
  
  cl_mem P_input; 
  cl_mem THETA_input;
  cl_mem IMG_output;
  
  size_t localworksize;
  size_t globalworksize;
  
  cl_event event;
  
  /*Kernel to execute */
  const char *source =   "\n" \
"__kernel void Backproject(const int N, const double ctr, const int numAngles,  \n" \
"                         const int projection_length, __global double *thetaPtr,  \n" \
"                         __global double *projections, __global double *img)  \n" \
"  {   \n" \
"  \n" \
"  int k, y;  \n" \
"  int input_row, output_row;  \n" \
"    \n" \
"  double xleft = -ctr;  \n" \
"  double ytop  = ctr;  \n" \
"  int center = projection_length/2;  \n" \
"  \n" \
"  int x = get_global_id(0);  \n" \
"  double xcoord = xleft + x;  \n" \
"  \n" \
"  double cos_theta, sin_theta;  \n" \
"  double t, fraction;  \n" \
"  int a;  \n" \
"  \n" \
"  if(x <= N)  \n" \
"  {  \n" \
"    output_row = x*N;  \n" \
"  \n" \
"    for(k=0;k<numAngles;++k)  \n" \
"    {  \n" \
"      cos_theta = cos(thetaPtr[k]);  \n" \
"      sin_theta = sin(thetaPtr[k]);  \n" \
"  \n" \
"      input_row = k * projection_length;  \n" \
"      t = xcoord * cos_theta + ytop * sin_theta;  \n" \
"  \n" \
"      for(y=0;y<N;++y)  \n" \
"      {  \n" \
"        a  = ((int) (t + N)) - N;  \n" \
"        fraction = t - a;  \n" \
"        a +=center;  \n" \
"  \n" \
"        img[output_row + y] += fraction * (projections[input_row + a + 1] - projections[input_row + a]) +  \n" \
"                               projections[input_row + a];  \n" \
"        t -= sin_theta;  \n" \
"      }  \n" \
"    }  \n" \
"  }  \n" \
"}  \n" \
"\n";
  
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
  
  projections = mxGetPr(P);
  projection_length = mxGetM(P);
  
  Nptr = mxGetPr(N_SIZE);
  N = (int) *Nptr;
  
  interp_ptr = mxGetPr(INTERP);
  interp_flag = (int) *interp_ptr;
  
  /* Create a matrix for the return argument */
  IMG = mxCreateDoubleMatrix(N, N, mxREAL);
  imgPtr = mxGetPr(IMG);
  
  center = (N-1) / 2;  
  xleft = -center;
  ytop = center;

  /* Initialization */
  error = clGetPlatformIDs(1, &platform, &no_platforms);
  error = clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 1, &device, &no_devices);
  error = clGetDeviceInfo(device, CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(size_t), &no_workgroups, NULL);
    
  context = clCreateContext(0, 1, &device, NULL, NULL, &error);
  commandqueue = clCreateCommandQueue(context, device, 0, &error);
  
  program_Backproject = clCreateProgramWithSource(context, 1, (const char **)&source, 
                                                    NULL, &error);
    
  error = clBuildProgram(program_Backproject, 0, NULL, NULL, NULL, NULL);
  
  kernel_Backproject = clCreateKernel(program_Backproject, "Backproject", &error);
  
  /* Localworksize dependent on if CPU or GPU is used. 1 for CPU, more give
   * no performance increase. N = 511 gives uneven number so round up to
   * 512 and handle extra in kernel. */
  localworksize = 16;
  globalworksize = N + 1;
  
  THETA_input  = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,  sizeof(double) * numAngles, thetaPtr, NULL);
  P_input      = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,  sizeof(double) * numAngles* projection_length, projections, NULL);
  IMG_output   = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(double) * N*N, NULL, NULL);
  error  = clSetKernelArg(kernel_Backproject, 0, sizeof(int), &N);
  error |= clSetKernelArg(kernel_Backproject, 1, sizeof(double), &center);
  error |= clSetKernelArg(kernel_Backproject, 2, sizeof(int), &numAngles);
  error |= clSetKernelArg(kernel_Backproject, 3, sizeof(int), &projection_length);
  error |= clSetKernelArg(kernel_Backproject, 4, sizeof(cl_mem), &THETA_input);
  error |= clSetKernelArg(kernel_Backproject, 5, sizeof(cl_mem), &P_input);
  error |= clSetKernelArg(kernel_Backproject, 6, sizeof(cl_mem), &IMG_output);

  /* Enqueue and run */
  error = clEnqueueNDRangeKernel(commandqueue, kernel_Backproject, 1, NULL, &globalworksize, &localworksize, 0, NULL, &event);
  clWaitForEvents(1, &event);
  
  /* Read result */
  error = clEnqueueReadBuffer(commandqueue, IMG_output, CL_TRUE, 0, sizeof(double) * N*N, imgPtr, 0, NULL, NULL); 
  clWaitForEvents(1, &event);
  
  clReleaseMemObject(IMG_output);
  clReleaseKernel(kernel_Backproject);
  clReleaseProgram(program_Backproject);
  clReleaseCommandQueue(commandqueue);
  clReleaseContext(context);
}
