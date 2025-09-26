function yf = makefilt(cut,len,ramp)

%   YF = MAKEFILT(CUT,LEN,RAMP) specifies a low-pass filter YF (a vector of
%   length LEN) in the spectral domain at the cutoff frequency CUT
%   with a transition ramp of width RAMP.
%
%   Input parameters:
%     CUT = cutoff frequency in terms of the sampling freq
%     LEN = length of the data to be filtered
%     RAMP = width of the transition ramp
% 
%   See also DATAFILT.
%

if cut >= 0.5
  error('Cutoff frequency too high.  Must be < 0.5');
end

if cut + ramp >= 0.5
  error('Ramp too wide.  Cut + ramp must be < 0.5');
end

high1=len/2;
high2=len/2-1;
if floor(len/2)~= len/2
  high2 = floor(len/2);
end

q = ( [ 0:high1 high2:-1:1 ]/len)';
if( ramp ~= 0 )
  t = (cut + ramp -q)/ ramp;
  t = t.*(t<0.99999).*(t>0.0001);
else
  t = 0*q;
end
yf = t + (q<=cut);

return


