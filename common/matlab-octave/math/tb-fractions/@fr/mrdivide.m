function FRAC = mrdivide(FRAC1,FRAC2)
% fr/mrdivide: matrix division for fraction objects
% usage: FRAC = FRAC1/FRAC2
%        FRAC = mrdivide(FRAC1,FRAC2)
% 
%
% arguments: (input)
%  FRAC1,FRAC2 - fraction arrays
%
% Note: right division by a matrix is currently not supported.
%
%
% arguments: (output)
%  FRAC - a fraction that represents FRAC1/FRAC2
%
%
% Example:
%  f = frac(1,2);
%  f/2
%  ans =
%     1 / 6
%
%
%  See also: rdivide, times, mtimes

% Author: Ben Petschel 25/7/09
%
% Version history:
%   25/7/09 - first release (using vpi/mrdivide as a template)

if (nargin~=2)
  error('fr:mrdivide:nargin','mrdivide must have two arguments');
end

% make sure both FRAC1 and FRAC2 are fraction objects
if ~isa(FRAC1,'fr')
  FRAC1 = fr(FRAC1);
end
if ~isa(FRAC2,'fr')
  FRAC2 = fr(FRAC2);
end

if (numel(FRAC2) == 1)
  % right division by a scalar
  FRAC = FRAC1./FRAC2;

else
  % no other modes are supported
  error('Sorry, matrix divide is not supported for this shape input arguments')
end

