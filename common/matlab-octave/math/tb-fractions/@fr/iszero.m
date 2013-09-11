function tf = iszero(FRAC)
% fr/iszero: test to see if a fraction object is zero
% usage: tf = iszero(FRAC);
% 
% arguments:
%  FRAC - a fraction object
%
%  tf  - returns a boolean, true or false.
%        tf == true if FRAC was zero
%
%
% Example:
%  iszero(fr(1,2))
%  ans = 
%     0
%
%  iszero(fr(0))
%  ans = 
%     1
%
%
%  See also: isunit

% Author: Ben Petschel 25/7/09
%
% Version history:
%   25/7/09 - first release (using vpi/iszero as a template)


if nargin == 0
  error('No input argument supplied')
end

if numel(FRAC)==1,
  tf = FRAC.whole==0 && FRAC.numer==0 && FRAC.denom==1;
  
else

  tf=false(size(FRAC));

  tf(:) = ([FRAC(:).whole]==0) & ([FRAC(:).numer]==0) & ([FRAC(:).denom]==1);

end;

end
