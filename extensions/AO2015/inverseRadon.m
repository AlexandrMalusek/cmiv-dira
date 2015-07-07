function [img,H] = inverseRadon(p, degVec, interpolation, filter, frequencyScaling, N)
	 
  global useCode;
  
  theta = pi*degVec/180;     %Convert degree to radian
  % Design the filter
  len=size(p,1);   
  H = designFilter(len, frequencyScaling);
  p(length(H),1)=0;  % Zero pad projections 
  
  p = fft(p);    % p holds fft of projections
  
  for i = 1:size(p,2)
    p(:,i) = p(:,i).*H; % frequency domain filtering
  end
  
  p = real(ifft(p));     % p is the filtered projections
  p(len+1:end,:) = [];   % Truncate the filtered projections
  %img = zeros(N,class(p));        % Allocate memory for the image.
  
  % Zero pad the projections to size 1+2*ceil(N/sqrt(2)) if this
  % quantity is greater than the length of the projections
  imgDiag = 2*ceil(N/sqrt(2))+1;  % largest distance through image.
  if size(p,1) < imgDiag 
    rz = imgDiag - size(p,1);  % how many rows of zeros
    p = [zeros(ceil(rz/2),size(p,2)); p; zeros(floor(rz/2),size(p,2))];
  end
  
  useCode = 2;
  
  % Backprojection - vectorized in (x,y), looping over theta
  switch useCode
    case 1
      img = Backprojectc(p, theta, N, 1);
    case 2
      img = Backprojectc_openmp(p, theta, N, 1);
    case 3
      img = Backprojectc_opencl(p, theta, N, 1);        
  end
  img = img*pi/(2*length(theta));
  
  function filt = designFilter(len, d)
	   
    order = max(64, 2^nextpow2(2*len));
    
    n = 0:(order/2); % 'order' is always even. 
    filtImpResp = zeros(1,(order/2)+1); % 'filtImpResp' is the bandlimited ramp's impulse response (values for even n are 0)
    filtImpResp(1) = 1/4; % Set the DC term 
    filtImpResp(2:2:end) = -1 ./ ((pi*n(2:2:end)).^2); % Set the values for odd n
    filtImpResp = [filtImpResp filtImpResp(end-1:-1:2)]; 
    filt = 2 * real(fft(filtImpResp)); 
    filt = filt(1:(order/2)+1);
    
    w = 2*pi*(0:size(filt, 2)-1)/order;   % frequency axis up to Nyquist
    
    filt(2:end) = filt(2:end) .* (1+cos(w(2:end)./d)) / 2;

    switch filter
      case 'cosine'
        filt(2:end) = filt(2:end) .* cos(w(2:end)/(2*d));
      case 'hann'
        filt(2:end) = filt(2:end) .* (1+cos(w(2:end)./d)) / 2;
      case 'hamming'
        filt(2:end) = filt(2:end) .* (.54 + .46 * cos(w(2:end)/d));
      case 'ram-lak'
         % No filtering
      case 'shepp-logan'
        filt(2:end) = filt(2:end) .* (sin(w(2:end)/(2*d)) ./ (w(2:end)/(2*d)));
      otherwise
        error(message('inverseRadon:invalidFilter'))
    end

    filt(w > pi*d) = 0;                      % Crop the frequency response
    filt = [filt' ; filt(end-1:-1:2)'];    % Symmetry of the filter
