function F=factor(FRAC)
% fr/factor - factors a fraction into a product of "prime" fractions
% usage: F=factor(FRAC)
%
% arguments:
%   FRAC - a scalar fraction object
%   F - a fraction object array where each element is a prime or the
%   inverse of a prime, and FRAC=prod(F)
%
% Example:
%  factor(fr(5,12))
%  ans =
%   1x4 fraction array with elements (reading down columns)
%     5
%     1 / 2
%     1 / 2
%     1 / 3
%
% see also: partial, prod

% Author: Ben Petschel 28/7/2009
%
% Version history:
%   28/7/2009 - first release

if nargin~=1
  error('fr:factor:nargin','factor only takes one argument');
end;
if numel(FRAC)~=1
  error('fr:factor:insize','factor only operates on scalars (non-arrays)');
end;
if ~isa(FRAC,'fr')
  try
    % convert to fraction
    FRAC=fr(FRAC);
  catch
    error('fr:factor:inclass','factor only operates on fraction-compatible objects');
  end;
end;

if ~isfinite(FRAC)
  % can't separate inf or nan
  F=FRAC;
elseif FRAC==0
  % leave 0 as-is
  F = FRAC;
  
else
  % find numerator and denominator of improper fraction and factor them
  [N,D]=rat(FRAC);
  if N==1
    fN = [];
  else
    fN = factor(N);
  end;
  if D==1
    fD = [];
  else
    fD = factor(D);
  end;
  F = [fr(0,fN,1),fr(0,1,fD)];
end;

