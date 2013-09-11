function FRAC = fix(FRAC)
% fr/fix: fix applied to a fraction object
%
% arguments:
%  FRAC  - a fraction object
%
%  See also: ceil, floor, round

% Author: Ben Petschel 25/7/09
%
% Version history:
%   25/7/09 - first release (using vpi/fix as a template)

% assume FRAC is reduced
for i=1:numel(FRAC)
  if FRAC(i).denom>0 && FRAC(i).numer>0
    % leave inf, nan, whole numbers unchanged
    FRAC(i).numer=0;
    FRAC(i).denom=1;
    K=FRAC(i).whole;
    if K<0
      % round numbers towards zero
      FRAC(i).whole=K+1;
    end;
  end;
end;
