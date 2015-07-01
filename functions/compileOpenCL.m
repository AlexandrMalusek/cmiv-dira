if(ispc)
  disp('Compiling OpenCL for Windows');
%  mex Backprojectc_opencl.cpp -I'C:\Program Files (x86)\AMD APP SDK\2.9-1\include' -L'C:\Program Files (x86)\AMD APP SDK\2.9-1\lib\x86_64' -lOpenCL
  mex computePolyProjc_opencl.cpp -I'C:\Program Files (x86)\AMD APP SDK\2.9-1\include' -L'C:\Program Files (x86)\AMD APP SDK\2.9-1\lib\x86_64' -lOpenCL
  mex sinogramJc_opencl.cpp -I'C:\Program Files (x86)\AMD APP SDK\2.9-1\include' -L'C:\Program Files (x86)\AMD APP SDK\2.9-1\lib\x86_64' -lOpenCL
elseif(isunix)
  disp('Compiling OpenCL for UNIX');
%  mex Backprojectc_opencl.cpp -lOpenCL -I'/opt/AMDAPPSDK-2.9-1/include' -L'/opt/AMDAPPSDK-2.9-1/lib/x86_64'
  mex computePolyProjc_opencl.cpp -lOpenCL -I'/opt/AMDAPPSDK-2.9-1/include' -L'/opt/AMDAPPSDK-2.9-1/lib/x86_64'
  mex sinogramJc_opencl.cpp -lOpenCL -I'/opt/AMDAPPSDK-2.9-1/include' -L'/opt/AMDAPPSDK-2.9-1/lib/x86_64'
end
