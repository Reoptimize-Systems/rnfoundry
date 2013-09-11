function B = le(FRAC1,FRAC2)
% fr/le: test for inequality (<=) between a pair of fraction objects
% usage: B = (FRAC1 <= FRAC2)
% usage: B = le(FRAC1,FRAC2)
% 
% arguments: (input)
%  FRAC1,FRAC2 - fraction objects or numeric values,
%                can be arrays provided the sizes are compatible
%
% arguments: (output)
%  result  - logical variable, true when FRAC1 <= FRAC2
%
%
%  See also: gt, lt, ge

% Author: Ben Petschel 25/7/09
%
% Version history:
%   25/7/09 - first release (using vpi/le as a template)


if (nargin ~= 2)
  error('a test for inequality must be between a pair of values')
end

B = sign(FRAC1-FRAC2)<=0;

end
