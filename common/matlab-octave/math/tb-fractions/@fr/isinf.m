function B=isinf(FRAC)
% fr/isinf: checks if a fraction object is infinity
%
% Usage: B=isinf(FRAC)
%
%  See also: isfinite, isnan

% Author: Ben Petschel 25/7/09
%
% Version history:
%   25/7/09 - first release

if nargin == 0
  error('No input argument supplied')
end

% infinite values have nonzero numerator and zero denominator
if numel(FRAC)==1,
  % use shortcuts in the scalar case
  B=(FRAC.denom==0)&&(FRAC.numer~=0);
else
  B=false(size(FRAC));
  B(:)=([FRAC(:).numer]~=0)&([FRAC(:).denom]==0);
end;
