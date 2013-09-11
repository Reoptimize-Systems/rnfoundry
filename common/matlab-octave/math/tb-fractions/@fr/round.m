function FRAC = round(FRAC)
% fr/round: round applied to a fraction object
%
% arguments:
%  FRAC  - a fraction object
%
%  See also: ceil, fix, floor

% Author: Ben Petschel 25/7/09
%
% Version history:
%   25/7/09 - first release (using vpi/round as a template)

for i=1:numel(FRAC)
  % assume FRAC is reduced
  N=FRAC(i).numer;
  D=FRAC(i).denom;
  if D>0 && N>0
    % leave inf, nan, whole numbers unchanged
    if N>=D-N
      % round up
      FRAC(i).whole=FRAC(i).whole+1;
    end;
    FRAC(i).numer=0;
    FRAC(i).denom=1;
  end;
end;