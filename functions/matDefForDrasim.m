function str = matDefForDrasim(matObj)
  % Return a string containing the definition of a material in Drasim
  %
  % matObj   Material object
  %
  % Example:
  % >> mat  = Material('adipose', 0.960, 'C71.3O16.4H10.0N1.8Ca0.3Cl0.2', 'maFr');
  % >> disp(matDefForDrasim(mat));
  % #define adipose formula=H0.582774C0.348695N0.007549O0.060211Cl0.000331Ca0.000440 rho=9.600000e-01

  % Compute normalized atomic fractions 
  n = length(matObj.W);
  atFr = zeros(n,1);
  sum = 0.0;
  for i = 1:n
    if matObj.W(i) ~= 0.0
      atFr(i) = matObj.W(i) / matObj.Ar(i);
      sum = sum + atFr(i);
    end
  end
  atFr = atFr / sum;   % Normalize the atomic fractions
  % Construct the string H0.582774C0.348695...
  atComp ='';
    for i = 1:n
    if matObj.W(i) ~= 0.0
       atComp = sprintf('%s%s%f', atComp, string(matObj.chemSymbol(i)),...
         atFr(i)); 
    end
  end

  str = sprintf('#define %s formula=%s rho=%e',...
    matObj.nameStr, atComp, matObj.density);
end
