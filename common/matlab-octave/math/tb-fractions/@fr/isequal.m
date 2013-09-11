function tf = isequal(FRAC1,FRAC2)
% fr/isequal: test to see if two fraction objects are identical
% usage: tf = isequal(FRAC1,FRAC2);
% 
% arguments:
% arguments: (input)
%  FRAC1,FRAC2 - fraction objects or compatible objects,
%                can be arrays of any size
%
% arguments: (output)
%  result  - logical variable, true when FRAC1 same as FRAC2
%
%
% Example:
%  isequal(fr(1,2),fr(1,3))
%  ans = 
%     0
%
%  isequal(fr(1),2)
%  ans = 
%     0
%
%  isequal(fr(1,7),1/7)
%  ans = 
%     1
%
%  isequal(fr(1,2),fr(0,1,2))
%  ans =
%     1
%
%
%  See also: iszero, isunit

% Author: Ben Petschel 25/7/09
%
% Version history:
%   25/7/09 - first release (using vpi/isequal as a template)


% a simple test.
tf = isequal(struct(fr(FRAC1)),struct(fr(FRAC2)));
