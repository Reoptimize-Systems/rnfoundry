function B=isnan(FRAC)
% fr/isnan: checks if a fraction object is not-a-number
%
% Usage: B=isnan(FRAC)
%
%  See also: isfinite, isinf

% Author: Ben Petschel 25/7/09
%
% Version history:
%   25/7/09 - first release

% non-a-numbers have zero numerator and denominator

if nargin == 0
  error('No input argument supplied')
end

if numel(FRAC)==1,
  % use shortcuts in the scalar case
  B=(FRAC.numer==0) && (FRAC.denom==0);
else
  B=false(size(FRAC));
  B(:)=([FRAC(:).numer]==0)&([FRAC(:).denom]==0);
end;
