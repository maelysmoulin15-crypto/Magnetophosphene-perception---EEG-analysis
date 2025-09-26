function y = datafilt(data,filt)

% Author: Alexandre Legros
% Copyright 2009 Alexandre Legros 

%   y = DATAFILT(DATA,FILT) filters the data in vector DATA in the
%   spectral domain with the filter described by the vector FILT
%   (created in MAKEFILT) to create the filtered data Y.
%
%   Note: the vectors DATA and FILT must have the same length.
%
%   See also MAKEFILT.
%
data = data(:);
if length(data) ~= length(filt)
  error( 'data and filt must have the same length.');
end
d = fft(data);
d = d.*filt;
y = real(ifft(d));

return

