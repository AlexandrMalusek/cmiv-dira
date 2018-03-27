/**
 * Computes polychromatic projections by summing the values of the 
 * projections of individual base materials over all photon energies used.
 * 
 * 
 */

#include "mex.h"
#include <math.h>

static void 
computePolychromaticProjection(int *ePtr, double ue, double *nPtr, double *pPtr,
                               double *muPtr, double *apPtr, int e_Size, int no_projections,
                                int mu_Size, int N, int M);

/* Input Arguments */
#define E     (prhs[0])
#define UE    (prhs[1])
#define N_P   (prhs[2])
#define P     (prhs[3])
#define MU    (prhs[4])

/* Output Arguments */
#define  AP    (plhs[0])

/**
 * Need 5 input arguments: E, uE, N, p, mu
 * Produces 1 output: Ap
 */
void 
mexFunction(int nlhs, mxArray  *plhs[], int nrhs, const mxArray  *prhs[])
{
  /* Input data */
  int *ePtr;                 /* Energies */
  double ue;  
  double *nPtr;              /* relative number of photons for Ul */
  double *pPtr;    
  double *muPtr;    
  
  /* temp pointers for copying values */
  double *from_doublePtr;   
  double *to_doublePtr;   
  int *from_intPtr;   
  int *to_intPtr;  
  int *dimPtr;   
  
  /*Size variables */
  int N;
  int M;
  int total_size;           /* Total image size, N * M */
  int e_Size;               /* Number of energies */
  int no_projections;       /* number of materials */
  int mu_Size;
  
  int k;                    /* Loop counter */

  /* Check validity of arguments */
  if (nrhs != 5)
  {
      mexErrMsgTxt("Incorrect number of INPUT arguments.");
  }
  
  if (nlhs != 1)
  {
      mexErrMsgTxt("Incorrect number of OUTPUT arguments.");
  }
  
  if (mxIsSparse(E) || mxIsSparse(UE) || mxIsSparse(N_P) || mxIsSparse(P) || mxIsSparse(MU))
  {
      mexErrMsgTxt("Sparse inputs not supported.");
  }
  
  if (!mxIsDouble(E) || !mxIsDouble(UE) || !mxIsDouble(N_P) || !mxIsDouble(P) || !mxIsDouble(MU))
  {
      mexErrMsgTxt("Input must be double.");
  }
  
  /* Matrix allocation */
  
  /* Get the size of E matrix and allocate memory */
  N = mxGetN(E);
  M = mxGetM(E);
  total_size = M*N; 
  ePtr = (int *) mxCalloc(total_size, sizeof(int));
  from_doublePtr = mxGetPr(E);
  to_intPtr = ePtr;
  for (k = 0; k < total_size; k++)
    *(to_intPtr++) = (int)*(from_doublePtr++);
  
  e_Size = M;
  
  /* Copy value of uE */
  from_doublePtr = mxGetPr(UE);
  ue = *from_doublePtr;
  
  /* Get the size of N matrix and allocate memory */
  N = mxGetN(N_P);
  M = mxGetM(N_P);
  total_size = M*N; 
  nPtr = (double *) mxCalloc(total_size, sizeof(double));
  from_doublePtr = mxGetPr(N_P);
  to_doublePtr = nPtr;
  for (k = 0; k < total_size; k++)
    *(to_doublePtr++) = *(from_doublePtr++);
  
  /* Get the size of MU matrix and allocate memory */
  N = mxGetN(MU);
  M = mxGetM(MU);
  total_size = M*N; 
  muPtr = (double *) mxCalloc(total_size, sizeof(double));
  from_doublePtr = mxGetPr(MU);
  to_doublePtr = muPtr;
  for (k = 0; k < total_size; k++)
    *(to_doublePtr++) = *(from_doublePtr++);
  
  mu_Size = M;
  
  /* P */
  N = mxGetN(P);
  M = mxGetM(P);
  total_size = M*N; 
  pPtr = (double *) mxCalloc(total_size, sizeof(double));
  from_doublePtr = mxGetPr(P);
  to_doublePtr = pPtr;
  for (k = 0; k < total_size; k++)
    *(to_doublePtr++) = *(from_doublePtr++);
  
  /* Get the z-dimension for matrix P */
  dimPtr = mxGetDimensions(P);
  no_projections = dimPtr[2];
  
  /* Allocate a 2D matrix for the output AP */
  /* Columns is 720x5 = 3600, need only 720 as column value, so divide by no_projections */
  AP = mxCreateDoubleMatrix(M, N/no_projections, mxREAL);
  
  computePolychromaticProjection(ePtr, ue, nPtr, pPtr, muPtr, mxGetPr(AP),
                                 e_Size, no_projections, mu_Size, N/no_projections, M);
}

static void 
computePolychromaticProjection(int *ePtr, double ue, double *nPtr, double *pPtr,
                               double *muPtr, double *apPtr, int e_Size, int no_projections,
                               int mu_Size, int N, int M)
{    
    /* Loop variables */
    int x,y;
    int k;
    int l;
    
    int energy;
    int image_size;
    
    double temporarySum;
    double result;
    
    int value;
    int smallest_x, smallest_y, largest_x, largest_y;
    
    image_size = M*N;
    
    /* Calculate for each pixel in the matrix */
    for(y=0;y<M;++y)
    {
        for(x=0;x<N;x++)
        {
            result = 0;

            for(k=0;k<e_Size;++k)
            {
                /* tmpSum = zeros(size(p(:, :, 1))); % 511x720 */
                temporarySum = 0;
                
                energy = ePtr[k];
                
                /* tmpSum = tmpSum+(-mu(E(k), i)*100.*p(:, :, i)); */
                for(l=0;l<no_projections;++l)
                {
                    temporarySum += -muPtr[l*mu_Size + energy - 1]*100*
                                     pPtr[y*N + x + l*image_size] ;
                }
                
                /* sl(:, :, k) = (E(k)*N(k)) .* exp(tmpSum);    */
                result += (energy * nPtr[k])*exp(temporarySum);
            }
            /* Ap = -log(up/uE);  */
            apPtr[y*N + x] = -log((result)/ue);
        }  
    }    
}   
