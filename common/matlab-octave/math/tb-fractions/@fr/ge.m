function B = ge(FRAC1,FRAC2)
% fr/ge: test for inequality (>=) between a pair of fraction objects
% usage: B = (FRAC1 >= FRAC2)
% usage: B = ge(FRAC1,FRAC2)
% 
% arguments: (input)
%  FRAC1,FRAC2 - fraction objects or compatible objects,
%                can be arrays provided the sizes are compatible
%
% arguments: (output)
%  result  - logical variable, true when FRAC1 >= FRAC2
% 
% Example:
%  FRAC1 = fr(1/2);
%  FRAC1 >= FRAC1
%  ans =
%     1
%
%  FRAC1 >= 1/3
%  ans =
%     1
%
%  FRAC1 >= 2/3
%  ans =
%     0
%
%
%  See also: gt, lt, le

% Author: Ben Petschel 25/7/09
%
% Version history:
%   25/7/09 - first release (using vpi/ge as a template)


if (nargin ~= 2)
  error('a test for inequality must be between a pair of values')
end

B = sign(FRAC1-FRAC2)>=0;

end


