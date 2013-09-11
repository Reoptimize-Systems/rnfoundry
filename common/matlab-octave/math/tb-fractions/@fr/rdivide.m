function FRAC = rdivide(FRAC1,FRAC2)
% fr/rdivide: quotient of two fraction objects
% usage: FRAC = FRAC1./FRAC2
% usage: FRAC = rdivide(FRAC1,FRAC2)
% 
%
% arguments: (input)
%  FRAC1,FRAC2 - fraction objects or compatible objects
%      can be arrays if the sizes are compatible
%
%
% arguments: (output)
%  FRAC - a fraction that represents the quotient FRAC1./FRAC2
%
%
% Example:
%  f = frac(1,2);
%  f/2
%  ans =
%     1 / 6
%
%
%  See also: mrdivide, times, mtimes

% Author: Ben Petschel 25/7/09
%
% Version history:
%   25/7/09 - first release (using vpi/rdivide as a template)

if (nargin~=2)
  error('fr:rdivide:nargin','rdivide must have two arguments');
end

% make sure both FRAC1 and FRAC2 are fraction objects
if ~isa(FRAC1,'fr')
  FRAC1 = fr(FRAC1);
end
if ~isa(FRAC2,'fr')
  FRAC2 = fr(FRAC2);
end

% multiply by inverse
FRAC = FRAC1.*frinv(FRAC2);

end
