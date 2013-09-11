function FRAC = floor(FRAC)
% fr/floor: floor applied to a fraction object
%
% arguments:
%  FRAC  - a fraction object
%
%  See also: ceil, fix, round

% Author: Ben Petschel 25/7/09
%
% Version history:
%   25/7/09 - first release (using vpi/floor as a template)

% assume FRAC is reduced
for i=1:numel(FRAC)
  if FRAC(i).denom>0 && FRAC(i).numer>0
    % leave inf, nan, whole numbers unchanged
    FRAC(i).numer=0;
    FRAC(i).denom=1;
  end;
end;
