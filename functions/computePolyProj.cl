#define energies 140 
#define images 5

kernel void computePolyProj(__global int *ePtr, const double ue, __global double *nPtr, 
                              __global double *pPtr, __global double *muPtr, 
                              __global double *apPtr, const int e_Size, 
                                const int no_projections, const int mu_Size, 
                                const int total_size) 
{ 
    int t_id = get_global_id(0); 
  
    int k; 
    int l; 
 
    __local int energyLevels[energies]; 
 
    for(k=1;k<energies;++k) 
        energyLevels[k] = ePtr[k]; 
   
    double temporarySum = 0; 
    double result = 0; 
 
    for(k = 1; k < e_Size-1; ++k) 
    { 
     temporarySum = 0; 
                
     for(l = 0; l < no_projections; ++l) 
      { 
       temporarySum += -muPtr[l*mu_Size + energyLevels[k] - 1]* 
                         100*pPtr[l*total_size + t_id]; 
      } 
     result += (energyLevels[k] * nPtr[k]) * 
                (energyLevels[k+1] - energyLevels[k-1]) * 
              exp(temporarySum); 
    } 
    apPtr[t_id] = -log(result/2/ue); 
} 