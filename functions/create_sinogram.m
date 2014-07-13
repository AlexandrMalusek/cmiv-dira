function [ im ] = create_sinogram(run, projections, dec_elements)
  %CREATE_SINOGRAM(run, projections, dec_elements)
  %
  %   Takes a set of projections created by drasim and creates a sinogram 
  %   of it.
  %
  %   run             = # of the run in drasim
  %   projections     = # of projections used in drasim
  %   dec_elements    = # of detector elements used in drasim
  %
  %   Oscar Grandell 2012
	 
  im = zeros(projections,dec_elements);
  fid = fopen(['bilder/b_d.',num2str(run),'.ima'], 'r');
  im(1,:) = fread(fid, [1,dec_elements], 'float32');
  fclose(fid);
  
  for i=1:projections-1
    fid = fopen(['bilder/b',num2str(i),'_d.',num2str(run),'.ima'], 'r');
    im(i+1,:) = fread(fid, [1,dec_elements], 'float32');
    fclose(fid);
  end
end
