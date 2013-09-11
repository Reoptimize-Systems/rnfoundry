function FRAC = abs(FRAC)
% fr/abs: absolute value of a fraction object
% 
% arguments:
%  FRAC - a fraction object
%
%  See also: sign

% Author: Ben Petschel 25/7/09
%
% Version history:
%   25/7/09 - first release (using vpi/abs as a template)

for i=1:numel(FRAC)
  K=FRAC(i).whole;
  if K<0
    FRAC(i).whole=-K;
    FRAC(i).numer=-FRAC(i).numer;
    FRAC(i)=freduce(FRAC(i));
  end;
end;
