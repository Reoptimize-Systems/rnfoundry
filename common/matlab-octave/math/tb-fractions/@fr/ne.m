function B = ne(FRAC1,FRAC2)
% fr/ne: test for inequality between a pair of fraction objects
% usage: B = (FRAC1 ~= FRAC2)
% usage: B = ne(FRAC1,FRAC2)
%
% arguments: (input)
%  FRAC1,FRAC2 - fraction objects or compatible objects,
%                can be arrays provided the sizes are compatible
%
% arguments: (output)
%  B  - logical variable, true when the two inputs
%            do NOT represent the same value (treating NaN's as equal).
%
%
%  See also: eq, ge, gt, lt, le

% Author: Ben Petschel 25/7/09
%
% Version history:
%   25/7/09 - first release (using vpi/ne as a template)


if (nargin ~= 2)
  error('a test for inequality must be between a pair of values')
end

B=~(FRAC1==FRAC2);

end



