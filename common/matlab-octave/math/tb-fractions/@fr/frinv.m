function FRAC=frinv(FRAC)
% fr/frinv - determines element-wise inverse of a fraction
% usage: FRAC = frinv(FRAC);
% 
% arguments:
%  FRAC - fraction object (scalar or array)
%
% See also: rdivide, mrdivide, times

% Author: Ben Petschel 25/7/09
%
% Version history:
%   25/7/09 - first release
%   15/12/09 - bug fix (handles non-doubles correctly)


if nargin~=1,
  error('fr:frinv:nargin','must have at least one input argument');
end;

K=FRAC.whole;
N=FRAC.numer;
D=FRAC.denom;

% 1/(k+n/d)=d/(kd+n)

if isa(N,'double'),
  FRAC.whole=0;
else
  FRAC.whole=N-N; % creates a zero of same type as K,N,D
end;
FRAC.numer=D;
FRAC.denom=K.*D+N; % vectorized

FRAC=freduce(FRAC); % reduce to lowest common terms
