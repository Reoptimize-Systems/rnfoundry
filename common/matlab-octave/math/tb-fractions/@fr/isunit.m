function tf = isunit(FRAC)
% fr/isunit: test to see if a fraction object is +1 or -1
% usage: tf = isunit(FRAC);
% 
% arguments:
%  FRAC - a fraction object or array
%
%  tf  - returns a boolean, true or false.
%        tf == true if FRAC was equal to 1
%
%
% Example:
%  isunit(fr(1,2))
%  ans = 
%     0
%
%  isunit(fr(1,1))
%  ans = 
%     1
%
%
%  See also: iszero

% Author: Ben Petschel 25/7/09
%
% Version history:
%   25/7/09 - first release (using vpi/isunit as a template)


if nargin == 0
  error('No input argument supplied')
end

if numel(FRAC)==1,
  tf = FRAC.whole==1 && FRAC.numer==0 && FRAC.denom==1;
  
else

  tf=false(size(FRAC));

  tf(:) = ([FRAC(:).whole]==1) & ([FRAC(:).numer]==0) & ([FRAC(:).denom]==1);

end;


