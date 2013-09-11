function FRAC = ceil(FRAC)
% fr/ceil: ceil applied to a fraction object
%
% arguments:
%  FRAC  - a fraction object
%
%  See also: fix, floor, round

% Author: Ben Petschel 25/7/09
%
% Version history:
%   25/7/09 - first release (using vpi/ceil as a template)

% assume FRAC is reduced
for i=1:numel(FRAC)
  if FRAC(i).denom>0 && FRAC(i).numer>0
    % leave inf, nan, whole numbers unchanged
    FRAC(i).numer=0;
    FRAC(i).denom=1;
    FRAC(i).whole=FRAC(i).whole+1;
  end;
end;
