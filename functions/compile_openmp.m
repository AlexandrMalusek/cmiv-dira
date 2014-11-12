if(ispc)
  disp('Compiling for Windows');
  mex openmp_Backprojectc.c COMPFLAGS="/openmp $COMPFLAGS"
  mex openmp_computePolyProjc.c COMPFLAGS="/openmp $COMPFLAGS"
  mex openmp_sinogramJc.c COMPFLAGS="/openmp $COMPFLAGS"
elseif(isunix)
  disp('Compiling for UNIX');
  mex openmp_Backprojectc.c CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
  mex openmp_computePolyProjc.c CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
  mex openmp_sinogramJc.c CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
end

mex Backprojectc.c
mex computePolyProjc.c
mex MD2c.c
mex MD3c.c
mex sinogramJc.c
mex WBHCc.c
