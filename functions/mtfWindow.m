function yarr = mtfWindow(xarr)

yarr = xarr * 0;

load -ascii MTF.dat
xaxis = MTF(:,1);
yaxis = MTF(:,2);
siz = size(xaxis);
xlength = siz(1);
siz = size(xarr);
xarrlength = siz(1)*siz(2);

for j = 1:xarrlength,
  xdata = xarr(j);
  k = 1;
  while (k<=xlength) & (xdata>xaxis(k))
    k=k+1;
  end;
  
  if k==1
    yarr(j) = yaxis(1);
  elseif k == xlength+1
    yarr(j) = 0;
  else
    weight1  = (xaxis(k)- xdata)      / (xaxis(k) - xaxis(k-1));
    weight2  = (xdata   - xaxis(k-1)) / (xaxis(k) - xaxis(k-1));
    yarr(j) = weight1 * yaxis(k-1) + weight2 * yaxis(k);
  end;
end;
