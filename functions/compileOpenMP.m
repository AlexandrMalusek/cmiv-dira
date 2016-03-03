if(ispc)
  disp('Compiling for Windows');
%  mex Backprojectc_openmp.c COMPFLAGS="/openmp $COMPFLAGS"
  mex computePolyProjc_openmp.c COMPFLAGS="/openmp $COMPFLAGS"
  mex sinogramJc_openmp.c COMPFLAGS="/openmp $COMPFLAGS"
elseif(isunix)
  disp('Compiling for UNIX');
%  mex Backprojectc_openmp.c CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
  mex computePolyProjc_openmp.c CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
  mex sinogramJc_openmp.c CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
end

%mex Backprojectc.c
mex computePolyProjc.c
mex MD2c.c
mex MD3c.c
mex sinogramJc.c
%mex WBHCc.c
