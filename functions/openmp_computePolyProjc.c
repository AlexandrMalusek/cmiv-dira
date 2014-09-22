#include "mex.h"
#include <math.h>
#include <omp.h>

static void 
computePolychromaticProjection(int *ePtr, double ue, double *nPtr, double *pPtr,
                               double *muPtr, double *apPtr, int e_Size, int p_Size,
			       int mu_Size, int columns, int rows);

static char rcs_id[] = "$Revision: 1.10 $";

/* Input Arguments */
#define E     (prhs[0])
#define UE    (prhs[1])
#define N     (prhs[2])
#define P     (prhs[3])
#define MU    (prhs[4])

/* Output Arguments */
#define	AP    (plhs[0])

/**
 * Need 5 input arguments: E, uE, N, p, mu
 * Produces 1 output: Ap
 */
void 
mexFunction(int nlhs, mxArray  *plhs[], int nrhs, const mxArray  *prhs[])
{
  /* Input pointers */
  int *ePtr;     /* Energy levels */
  double ue;	
  double *nPtr;		
  double *pPtr;		
  double *muPtr;    
  
  /* temp pointers for copying values */
  double *from_doublePtr;   
  double *to_doublePtr;   
  int *from_intPtr;   
  int *to_intPtr;  
  
  /*Size variables */
  int columns;      /* columns */
  int rows;         /* rows */
  int total_size;   /* Total image size, columns * rows */
  int e_Size;       /* Number of energy levels */
  int p_Size = 5;   /* number of layers */
  int mu_Size;
  
  int k;            /* Loop counter */
  
  /* Check validity of arguments */
  if (nrhs != 5)
    {
      mexErrMsgTxt("Incorrect number of INPUT arguments.");
    }
  
  if (nlhs != 1)
    {
      mexErrMsgTxt("Incorrect number of OUTPUT arguments.");
    }
  
  if (mxIsSparse(E) || mxIsSparse(UE) || mxIsSparse(N) || mxIsSparse(P) || mxIsSparse(MU))
    {
      mexErrMsgTxt("Sparse inputs not supported.");
    }
  
  if (!mxIsDouble(E) || !mxIsDouble(UE) || !mxIsDouble(N) || !mxIsDouble(P) || !mxIsDouble(MU))
    {
      mexErrMsgTxt("Input must be double.");
    }
  
  /**
   * Matrix allocation
   */
  
  /* Get the size of E matrix and allocate memory */
  columns = mxGetN(E);
  rows = mxGetM(E);
  total_size = rows * columns; 
  ePtr = (int *) mxCalloc(total_size, sizeof(int));
  from_doublePtr = mxGetPr(E);
  to_intPtr = ePtr;
  for (k = 0; k < total_size; k++)
    *(to_intPtr++) = (int)*(from_doublePtr++);
  
  e_Size = rows;
  
  /* Copy value of uE */
  from_doublePtr = mxGetPr(UE);
  ue = *from_doublePtr;
  
  /* Get the size of N matrix and allocate memory */
  columns = mxGetN(N);
  rows = mxGetM(N);
  total_size = rows * columns; 
  nPtr = (double *) mxCalloc(total_size, sizeof(double));
  from_doublePtr = mxGetPr(N);
  to_doublePtr = nPtr;
  for (k = 0; k < total_size; k++)
    *(to_doublePtr++) = *(from_doublePtr++);
  
  /* Get the size of MU matrix and allocate memory */
  columns = mxGetN(MU);
  rows = mxGetM(MU);
  total_size = rows * columns; 
  muPtr = (double *) mxCalloc(total_size, sizeof(double));
  from_doublePtr = mxGetPr(MU);
  to_doublePtr = muPtr;
  for (k = 0; k < total_size; k++)
    *(to_doublePtr++) = *(from_doublePtr++);
  
  mu_Size = rows;
  
  /* P */
  columns = mxGetN(P);
  rows = mxGetM(P);
  total_size = rows * columns; 
  pPtr = (double *) mxCalloc(total_size, sizeof(double));
  from_doublePtr = mxGetPr(P);
  to_doublePtr = pPtr;
  for (k = 0; k < total_size; k++)
    *(to_doublePtr++) = *(from_doublePtr++);
  
  /* Allocate a 2D matrix for the output AP */
  AP = mxCreateDoubleMatrix(511, 720, mxREAL);
  
  /* Columns is 720x5 = 3600, need only 720 as column value, so divide by p_Size */
  computePolychromaticProjection(ePtr, ue, nPtr, pPtr, muPtr, mxGetPr(AP),
                                 e_Size, p_Size, mu_Size, columns/p_Size, rows);
}

static void 
computePolychromaticProjection(int *ePtr, double ue, double *nPtr, double *pPtr,
                               double *muPtr, double *apPtr, int e_Size, int p_Size,
			       int mu_Size, int columns, int rows)
{    
  /* Loop variables */
  int i;
  int j;
  int k;
  int l;
  
  int energyLevel;
  
  double temporarySum;
  double result;
  
#pragma omp parallel for private(j, k, l, energyLevel, temporarySum, result)
  for(i = 0; i < rows; ++i)
    {
      for(j = 0; j < columns; j++)
        {
	  result = 0;
	  
	  /* tmpSum = tmpSum+(-mu(E(k), i)*100.*p(:, :, i)); */
	  for(k = 1; k < e_Size-1; ++k)
            {
	      /* tmpSum = zeros(size(p(:, :, 1))); % 511x720 */
	      temporarySum = 0;
              
	      energyLevel = *(ePtr + k);
	      
	      for(l = 0; l < p_Size; ++l)
                {
		  temporarySum += (-*(muPtr + l*mu_Size + energyLevel - 1))*
		    100*(*(pPtr + i*columns + j + l*rows*columns));
                }
	      
	      /* sl(:, :, k) = (E(k)*N(k))*(E(k+1)-E(k-1)).*exp(tmpSum); */
	      result += (energyLevel * (*(nPtr + k)))*
		( *(ePtr + k + 1) - *(ePtr + k - 1))*
		exp(temporarySum);
            }
	  /* Ap = -log(up/uE); */
	  *(apPtr + i*columns + j) = -log((result/2)/ue);
        }  
    }
}
