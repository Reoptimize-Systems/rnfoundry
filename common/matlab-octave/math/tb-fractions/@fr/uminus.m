function FRAC = uminus(FRAC)
% fr/uminus: Negates a fraction object
% usage: FRAC = -FRAC;
%        FRAC = uminus(FRAC);
% 
% arguments:
%  FRAC - a fraction object
%
%
%  See also: minus, uplus

% Author: Ben Petschel 25/7/09
%
% Version history:
%   25/7/09 - first release (using vpi/uminus as a template)

for i=1:numel(FRAC)
  FRAC(i).whole = -FRAC(i).whole;
  FRAC(i).numer = -FRAC(i).numer;
end;

FRAC = freduce(FRAC);
