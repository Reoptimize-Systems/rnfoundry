function [V,rep]=cfracsqrt(f,a,b,m)
% fr/cfracsqrt: find the continued fraction of a*sqrt(f)+b
% usage: [V,rep]=cfracsqrt(f)
%        [V,rep]=cfracsqrt(f,a)
%        [V,rep]=cfracsqrt(f,a,b)
%        [V,rep]=cfracsqrt(f,a,b,m)
%
% CFRACSQRT returns the continued fraction expansion of a*sqrt(f)+b
% (a=1, b=0 by default).  If f is square, a*sqrt(f)+b is a fraction and the
% continued fraction terminates; otherwise the continued fraction will
% repeat after a few terms.
%
% CFRACSQRT determines the terms of the continued fraction until either the
% repeating part is found, or the denominator of the convergent is larger
% than m.  If the repeating part is found, the full continued fraction
% is infinite and of the form:
%   [V,rep,rep,...]
% If the denominator limit is reached first, then rep=[] is returned and
% V matches the first few terms of the continued fraction.
%
% Example: determine the best rational approximations of sqrt(13) with
% denominator at most 20:
%   f=fr(13);a=1;b=0;
%   dmax=20;
%   [cf,rep]=cfracsqrt(f)
%   [r1,r2]=bestrat(fr(cf),rep,dmax)
%   if sign(a)*(r1+r2)^2 < sign(a)*4*a^2*f, rbest=r2; else, rbest=r1; end;
%
%
% Theory:
%
% Continued fractions are defined recursively: if x = [n0,n1,...] where nk
% are all integers, then x = n0 + 1/x1 where x1 = [n1,...] > 1 (put []=inf
% and [nk]=nk).  So n0 = floor(x) and x1 = 1/(x-n0) if x>n0.
%
% At each step xk can be expressed as ak*sqrt(f)+bk where ak and bk are
% fractions and eventually [ak,bk] repeat.
% 
% The tricky part is calculating n=floor(a*sqrt(f)+b):
% effectively we maximize n such that n <= (a*sqrt(f)+b)
% This can be rearranged as a piecewise quadratic inequality in n and
% solved exactly (within the limits of the fraction components data types).
%
% See also: cfrac, bestrat


% Author: Ben Petschel 18/8/2009
% Version history:
%   18/8/2009 - first release



% at each stage we have term a*sqrt(f)+b
% put n=floor(a*sqrt(f)+b)
% -> get n+(a*sqrt(f)+b-n) = n+1/(1/(a*sqrt(f)+b-n))
%    = n+1/(a*sqrt(f)+n-b)/(a^2*f-(b-n)^2)
%    = n+1/((a/d)*sqrt(f)+(n-b)/d) where d=a^2*f-(b-n)^2

if nargin<1 || nargin>4
  error('fr:cfracsqrt:nargin','requires 1 to 4 arguments');
end;
if numel(f)>1
  error('fr:cfracsqrt:numel','f must be a scalar fraction');
end;
if ~isfinite(f)
  error('fr:cfracsqrt:finite','f must be finite');
end;
if f<0
  error('fr:cfracsqrt:negative','f must be non-negative');
end;
if nargin<2
  a=1;
end;
if nargin<3
  b=0;
end;
limdenom=(nargin==4);

if ~isa(f,'fr')
  f=fr(f);
end;
if ~isa(a,'fr')
  a=a+(f-f); % forces a to be same type as f
end;
if ~isa(b,'fr')
  b=b+(f-f); % forces b to be same type as f
end;

alist=a;
blist=b;

if limdenom
  % keep track of convergents
  h1=0;
  k1=1;
  h=1;
  k=0;
end;

keepgoing=true;
V=[];
while keepgoing
  n=floorsqrt(f,a,b);
  % check if denominator limit reached or repeating part has been found
  V=[V,n];
  if limdenom
    % check denominator
    h2=h1;
    k2=k1;
    h1=h;
    k1=k;
    h=h2+h1*n;
    k=k2+k1*n;
    if k>m
      keepgoing=false;
      rep=[]; % repeating part not found
    end;
  end;
  d=(a^2*f-(b-n)^2);
  if d==0
    % trying to calculate 1/0 (continued fraction terminates here)
    keepgoing = false;
    rep=[]; % no repeating part
  else
    a=a/d;
    b=(n-b)/d;
    ind=find(alist==a & blist==b,1,'first');
    if isempty(ind)
      alist=[alist,a];
      blist=[blist,b];
    else
      % found repeating part
      keepgoing=false;
      rep = V(ind:end);
      V = V(1:ind-1);
    end;
  end;
end;

end

function n=floorsqrt(f,a,b)
% finds n=floor(a*sqrt(f)+b)
% n <= (a*sqrt(f)+b) -> (n-b) <= a*sqrt(f)
% (i) a>0 -> 0 <= max(0,n-b) <= a*sqrt(f)
%         -> max(0,n-b)^2 <= a^2*f  (*)
% (ii) a<0 -> min(0,n-b) <= a*sqrt(f) <= 0
%         -> -min(0,n-b) >= -a*sqrt(f) >= 0
%         -> max(0,n-b) >= -a*sqrt(f) >= 0
%         -> max(0,b-n)^2 >= a^2*f  (**)
%
% Solve this by finding an interval [n1,n2] such that h(n1) and h(n2) are
% on opposite sides of a^2*f where h(n) is LHS of (*) or (**).
% Once the interval is found, bisect until n2-n1=1 or a solution is found
switch sign(a)
  case 0
    n=floor(b);
  case 1
    % a>0: maximize n s.t. max(0,n-b)^2 <= a^2*f
    % start with interval floor(double(a)*sqrt(double(f)+double(b))+[0,1]
    % in most cases this will contain (a*sqrt(f)+b) but occasionally not,
    % e.g. for very large numbers or vpi's
    try
      n1=(b-b)+floor(double(a)*sqrt(double(f))+double(b)); % forces n1 to fraction
    catch
      % above can fail when coefficient are vpi and n1>2^53-1
      n1=floor(b);
    end;
    n2=n1+1;
    rhs=(a^2)*f;
    if n1>b,n1rem=(n1-b)^2;else n1rem=0;end;
    %n1rem = max(0,n1-b)^2;
    if n2>b,n2rem=(n2-b)^2;else n2rem=0;end;
    %n2rem = max(0,n2-b)^2;
    while n1rem >= rhs
      % double interval (moving left) until found
      n2old=n2;
      n2=n1;
      n2rem=n1rem;
      n1=3*n1-2*n2old; % n1new = n2new+2*(n1old-n2old) = 3*n1-2*n2old
      if n1>b,n1rem=(n1-b)^2;else n1rem=0;end;
      %n1rem=max(0,n1-b)^2;
    end;
    while n2rem < rhs
      % double interval (moving right) until found
      n1old=n1;
      n1=n2;
      n1rem=n2rem;
      n2=3*n2-2*n1old; % n2new = n1new+2*(n2old-n1old) = 3*n2-2*n1old
      if n2>b,n2rem=(n2-b)^2;else n2rem=0;end;
      %n2rem=max(0,n2-b)^2;
    end;
    while (n2-n1>1) && (n2rem~=rhs)
      % bisect interval until n found (maximizing n s.t. max(0,n-b)^2<=rhs)
      nmid = floor((n1+n2)/2);
      if nmid>b,nmidrem=(nmid-b)^2;else nmidrem=0;end;
      %nmidrem = max(0,nmid-b)^2;
      if nmidrem >= rhs
        n2=nmid;
        n2rem=nmidrem;
      else
        n1=nmid;
        n1rem=nmidrem;
      end;
    end;
    if n2rem==rhs
      n=n2;
    else
      n=n1;
    end;
  otherwise
    % case a<0
    % maximize n s.t. max(0,b-n)^2 >= a^2*f
    % start with interval floor(double(a)*sqrt(double(f))+double(b))+[0,1]
    % on most occasions this will include (a*sqrt(f)+b) but sometimes not,
    % e.g. with very large numbers or vpi
    try
      n1=(b-b)+floor(double(a)*sqrt(double(f))+double(b)); % forces n1 to fraction
    catch
      % above can fail when coefficient are vpi and n1>2^53-1
      n1=floor(b);
    end;
    n2=n1+1; % have b-n2<0
    rhs=(a^2)*f;
    %n1rem = (b-n1)^2;
    %n2rem = 0;
    if n1<b,n1rem=(b-n1)^2;else n1rem=0;end;
    %n1rem = max(0,b-n1)^2;
    if n2<b,n2rem=(b-n2)^2;else n2rem=0;end;
    %n2rem = max(0,b-n2)^2;
    while n2rem >= rhs
      % double interval (moving right) until found
      n1old=n1;
      n1=n2;
      n1rem=n2rem;
      n2=3*n2-2*n1old; % n2new = n1new + 2*(n2old-n1old) = 3*n2-2*n1old
      if n2<b,n2rem=(b-n2)^2;else n2rem=0;end;
      %n2rem = max(0,b-n2)^2;
    end;
    while n1rem < rhs
      % double interval (moving left) until found
      n2old=n2;
      n2=n1;
      n2rem=n1rem;
      n1=3*n1-2*n2old; % n1new = n2new+2*(n1old-n2old) = 3*n1-2*n2old
      if n1<b,n1rem=(b-n1)^2;else n1rem=0;end;
      %n1rem=max(0,b-n1)^2;
    end;
    % now have interval [n1,n2] s.t. (b-n1)^2 >= rhs and (b-n2)^2 < rhs
    while (n2-n1>1) && (n1rem~=rhs)
      % bisect interval until n found (maximizing n s.t. max(0,b-n)^2>=rhs)
      nmid = floor((n1+n2)/2);
      if nmid<b,nmidrem=(b-nmid)^2;else nmidrem=0;end;
      %nmidrem = max(0,b-nmid)^2;
      if nmidrem >= rhs
        n1=nmid;
        n1rem=nmidrem;
      else
        n2=nmid;
        n2rem=nmidrem;
      end;
    end;
    if n1rem==rhs
      n=n1;
    else
      n=n2;
    end;
end; % switch statement
n=n.whole;

end % helper function floorsqrt(...)