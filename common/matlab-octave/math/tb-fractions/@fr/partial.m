function FPART = partial(FRAC,flag)
% fr/partial - separates FRAC into partial fractions
% usage: FPART = partial(FRAC)
%
% arguments:
%  FRAC - a scalar fraction object
%  flag - boolean indicating if the "full" partial expansion is required
%  FPART - an array of fraction objects of the form a/p^k
%          where the p are prime such that FPART = sum(FPART)
%          If flag = true or is not provided, finds 0<a<p (full expansion)
%          If flag = false, finds 0<a<p^k.  See examples below
%
%
% Note:
%  FRAC can be a fraction of any compatible objects, provided FACTOR and
%  GCD (3 output form) are defined.
%
% Example:
%  f = fr(1,12);  % 1 / 12
%  fp = partial(f,true)
%  fp =
%   1x4 fraction array with elements (reading down columns)
%     -1
%     1 / 2
%     1 / 4
%     1 / 3
%
%  sum(fp)
%  ans =
%     1 / 12
%
%  partial(f,false) % don't split terms a/p^k
%  ans =
%   1x3 fraction array with elements (reading down columns)
%     -1
%     3 / 4
%     1 / 3
%
%
% See also: fr

% Author: Ben Petschel 28/7/09
%
% Version history:
%   28/7/09 - first release
%   30/7/09 - added option "flag" for non-full expansion


if nargin<1 || nargin>2,
  error('fr:partial:nargin','partial only takes one or two arguments');
end;
if nargin==1
  flag = true;
else
  % convert flag to boolean if necessary
  if ~isscalar(flag)
    error('fr:partial:flagsize','flag must be a scalar variable');
  end;
  if ~islogical(flag)
    try
      flag = (flag~=0);
    catch
      error('fr:partial:flagtype','flag must be convertable to a logical');
    end;
  end;
end;
if numel(FRAC)~=1
  error('fr:partial:insize','partial only operates on scalars (non-arrays)');
end;
if ~isa(FRAC,'fr')
  try
    % convert to fraction
    FRAC=fr(FRAC);
  catch
    error('fr:partial:inclass','partial only operates on fraction-compatible objects');
  end;
end;

if ~isfinite(FRAC)
  % can't separate inf or nan
  FPART=FRAC;
elseif FRAC==0
  % leave 0 as-is
  FPART = FRAC;
  
else
  
  % determine the whole term and remainder
  FPART = floor(FRAC);
  FRAC = FRAC - FPART;
  
  if FRAC~=0
    f=factor(FRAC.denom);
  end;
  while ~isempty(f)
    % denominator is always > 1
    p=f(1);
    feqp=(f==p); % boolean index of factors matching f(1)
    k=sum(feqp); % number of times f(1) appears in denominator
    % split fraction into a/p^k + b/q
    if k==numel(f)
      % denominator is p^k; this is the last term to split
      Fpk = FRAC;
      FRAC = 0; % zero left over
      f=[];
    else
      % split p^k from remainder
      D1=p^k;
      D2=FRAC.denom/D1;
      [g,c,d]=gcd(D1,D2); % c*D1+d*D2=1, so d/D1+c/D2=1/(D1*D2)
      if g~=1
        % don't expect this error to occur:
        error('factor returned non-coprime factors');
      end;
      N = FRAC.numer;
      N1 = d*N;
      N2 = c*N;
      Fpk = fr(0,N1,D1);
      FRAC = fr(0,N2,D2); % could use FRAC-Fpk but this should be the same
      f(feqp)=[]; % removes p from factor list
    end;
    % move the whole parts to FPART(1), then split 
    Fpkwhole = floor(Fpk);
    FRACwhole = floor(FRAC);
    FPART(1)=FPART(1)+Fpkwhole+FRACwhole;
    Fpk=Fpk-Fpkwhole;
    FRAC=FRAC-FRACwhole;
    if flag
      % split Fpk=a/p^k into sum of b/p^j, 1<=j<=k
      dig=digits(Fpk,k,p); % digits gives expansion of Fpk into sum(a/p^j)
      FPART = [FPART,fr(dig,p.^(1:k))];
    else
      % keep Fpk as-is
      FPART = [FPART,Fpk];
    end;
  end;
  if numel(FPART)>1 && FPART(1)==0
    % remove zero term
    FPART = FPART(FPART~=0);
  end;
end;

