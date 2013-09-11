function B = lt(FRAC1,FRAC2)
% fr/lt: test for strict inequality (<) between a pair of fraction objects
% usage: B = (FRAC1 < FRAC2)
% usage: B = lt(FRAC1,FRAC2)
% 
% arguments: (input)
%  FRAC1,FRAC2 - fraction objects or numeric values,
%                can be arrays provided the sizes are compatible
%
% arguments: (output)
%  result  - logical variable, true when FRAC1 < FRAC2 (vectorized)
%
%
%  See also: gt, le, ge

% Author: Ben Petschel 25/7/09
%
% Version history:
%   25/7/09 - first release (using vpi/lt as a template)


if (nargin ~= 2)
  error('a test for inequality must be between a pair of values')
end

B = sign(FRAC1-FRAC2)<0;

end





