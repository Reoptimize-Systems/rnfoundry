function [V,C]=cfrac(F)
% fr/cfrac - determines the continued fraction of the fraction object F
% usage: V=cfrac(F)
%        [V,C]=cfrac(F)
% 
% CFRAC returns the coefficients of the continued fraction expansion of F.
% [V,C]=cfrac(F) also returns the convergents of the continued fraction.
%
% The convergents of a continued fraction are a subset of the best rational
% approximations of a number.  For the best rational approximation with a
% given denominator bound, use BESTRAT.
%
% Example: determine the convergents of 3/8:
%   [V,C]=cfrac(fr(3,8)) % returns [0,2,1,2] and [0,1/2,1/3,3/8]
%
% Example: evaluate continued fraction represented by [1,2,3,4]
%   bestrat(fr([1,2,3,4])) % returns 1+13/20
%
%
% Theory:
%
% Continued fractions are defined recursively: if x = [n0,n1,...] where nk
% are all integers, then x = n0 + 1/x1 where x1 = [n1,...] > 1 (put []=inf
% and [nk]=nk).  So n0 = floor(x) and x1 = 1/(x-n0) if x>n0.
%
%
% see also: bestrat, cfracsqrt

% Author: Ben Petschel 18/8/2009
% Version history:
%   18/8/2009 - first release

getconv=(nargout>1); % whether or not to get the convergents

V=F.whole;
F.whole = 0;

if getconv
  n2=0;
  d2=1;
  n1=1;
  d1=0;
  n=n2+n1*V;
  d=d2+d1*V;
  N=n; % vector of numerators in convergents
  D=d; % vector of denominators in convergents
end;

while F>0
  F=1/F;
  V=[V,F.whole];
  F.whole=0;
  if getconv
    n2=n1;
    d2=d1;
    n1=n;
    d1=d;
    n=n2+n1*V(end);
    d=d2+d1*V(end);
    N=[N,n];
    D=[D,d];
  end;
end;

if getconv
  C=fr(N,D);
end;
