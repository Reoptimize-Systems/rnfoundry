function [r1,r2]=bestrat(f,varargin)
% fr/bestrat: determines the best rational approximations to a number
% usage:       r=bestrat(cf)
%        [r1,r2]=bestrat(f,d)
%        [r1,r2]=bestrat(cf,d)
%        [r1,r2]=bestrat(cf,rep,d)
%
% bestrat(f,d) returns the best rational approximations to a fraction f
% such that r1<=f<=r2 and the denominators of r1 and r2 are at most d
%
% bestrat(cf,d), where cf is an array of fractional integers, returns the
% best rational approximation of the fraction represented by cf,
% i.e. f=cf(1)+1/(cf(2)+1/...)
%
% bestrat(cf,rep,d) returns the best rational approximation to the
% number with the infinite repeating continued fraction [cf,rep,rep,...].
%
% bestrat(cf) or bestrat(cf,inf) returns the fraction represented by cf
%
% Note: the first argument must be a fraction object; e.g. pass outputs of
% cfrac via "fr":  bestrat(fr(cf))
%
% The choice of which of r1 and r2 is the better rational approximation
% depends on the form of the original fraction.
% If cf represents f, then (r2-f) < (f-r1) if  (r1+r2) < 2*f
% If cf represents x=a*sqrt(f)+b, then (r2-x) < (x-r1) if (r1+r2) < 2*x
% i.e. (r1+r2-2*b) < 2*a*sqrt(f)
% i.e. sign(a)*(r1+r2-2*b) < sign(a)*4*a^2*f
% 
%
% Example: evaluate continued fraction represented by [1,2,3,4]
%  bestrat(fr([1,2,3,4])) % returns 1+13/20
%
% Example: determine the best rational approximations of sqrt(13) with
% denominator at most 20 (3+3/5):
%   f=fr(13);a=1;b=0;
%   dmax=20;
%   [cf,rep]=cfracsqrt(f)
%   [r1,r2]=bestrat(fr(cf),rep,dmax)
%   if sign(a)*(r1+r2-2*b)^2<sign(a)*4*a^2*f,rbest=r2,else,rbest=r1,end;
%
%
% see also cfrac, cfracsqrt

% Author: Ben Petschel 18/8/2009
% Version history:
%   18/8/2009 - first release

if numel(f)==1 && nargin==2
  % no repeating part was specified
  cf=cfrac(f);
  rep = [];
  dmax=varargin{1};
else
  cf=[f(:).whole];
  if any(cf-f~=0)
    error('fr:bestrat:wholefrac','bestrat currently accepts only fractions representing whole numbers');
  end;
  if nargin == 1
    rep = [];
    dmax = inf;
  elseif nargin == 2
    rep = [];
    dmax = varargin{1};
  elseif nargin == 3
    rep = varargin{1};
    if isa(rep,'fr'),
      if any(rep~=floor(rep)),
        error('fr:bestrat:wholefrac','bestrat currently accepts only fractions representing whole numbers');
      end;
      rep=[rep(:).whole];
    end;
    dmax = varargin{2};
  else
    error('fr:bestrat:nargin','bestrat takes at most 3 inputs');
  end;
end;
rfull = ~isa(dmax,'vpi')&&isinf(dmax); % if true, evaluate full cf

k=1;
n2=0;
d2=1;
n1=1;
d1=0;
n=n2+n1*cf(k);
d=d2+d1*cf(k); % initially have d=1

keepgoing = true;
if ~rfull && d>=dmax
  keepgoing=false;
  r1=fr(n,d);
  r2=fr(n1,d1);
end;

ncf = numel(cf);
nrep = numel(rep);

while keepgoing
  % add terms to continued fraction
  k=k+1;
  n2=n1;
  d2=d1;
  n1=n;
  d1=d;
  if k>ncf && nrep==0
    % exact approximation is r1=r2=f
    keepgoing=false;
    r1=fr(n1,d1);
    r2=r1;
  else
    if k>ncf
      cfk = rep(mod(k-ncf-1,nrep)+1);
    else
      cfk = cf(k);
    end;
    d=d2+d1*cfk;
    if rfull || d<dmax
      n=n2+n1*cfk;
    else
      % d=d2+d1*b where b>=cf(k)
      keepgoing=false;
      b=floor((dmax-d2)/d1);
      % interval is (n1/d1,n/d) or (n/d,n1/d1)
      d=d2+d1*b;
      n=n2+n1*b;
      if b==cfk && k==ncf && nrep==0
        %full cf has denominator <= dmax
        r1=fr(n,d);
        r2=r1;
      elseif mod(k,2)==0
        r1=fr(n1,d1);
        r2=fr(n,d);
      else
        r1=fr(n,d);
        r2=fr(n1,d1);
      end;
    end;
  end;
  %[vpi(n),vpi(d)]
end;

end
