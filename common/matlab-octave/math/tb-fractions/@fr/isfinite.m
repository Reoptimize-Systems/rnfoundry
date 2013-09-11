function B=isfinite(FRAC)
% fr/isfinite: checks if a fraction object is finite
%
% Usage: B=isfinite(FRAC)
%
%  See also: isnan, isinf

% Author: Ben Petschel 25/7/09
%
% Version history:
%   25/7/09 - first release

if nargin == 0
  error('No input argument supplied')
end

% finite numbers have nonzero denominator
if numel(FRAC)==1,
  % simpler referencing for scalar case
  B=FRAC.denom~=0;
else
  B=false(size(FRAC));
  B(:)=[FRAC(:).denom]~=0;
end;
