function FRAC = plus(FRAC1,FRAC2)
% fr/plus: Adds two fraction objects, or adds a numeric var to a fraction
% usage: FRAC = FRAC1 + FRAC2;
%
% arguments: (input)
%  FRAC1,FRAC2 - fraction objects or any compatible objects (e.g. double)
%                either can be arrays provided the sizes are compatible.
%
% arguments: (output)
%  FRAC - a fraction object representing the sum:
%       FRAC = FRAC1 + FRAC2
%
% Example:
%  f = fr(1,3);
%  f+1
%  ans =
%     1 + 1 / 3
%
%  fr(1,3)+fr(2,3)
%  ans = 
%     1
%
%
%  See also: minus, uplus

% Author: Ben Petschel 25/7/09
%
% Version history:
%   25/7/09 - first release (using vpi/plus as a template)
%   22/8/09 - slight speedup using mygcd (less argcheck overhead)

if nargin~=2
  error('Plus is a dyadic operator, exactly 2 arguments')
end

% make sure both FRAC1 and FRAC2 are fraction objects
if ~isa(FRAC1,'fr')
  FRAC1 = fr(FRAC1);
end
if ~isa(FRAC2,'fr')
  FRAC2 = fr(FRAC2);
end

n1 = numel(FRAC1);
n2 = numel(FRAC2);

if (n1 == 1) && (n2 == 1)
  % a pair of scalars

  % do the addition
  if isfinite(FRAC1) && isfinite(FRAC2)
    % adding two real numbers
    
    K1=FRAC1.whole;
    K2=FRAC2.whole;
    N1=FRAC1.numer;
    N2=FRAC2.numer;
    D1=FRAC1.denom;
    D2=FRAC2.denom;

    % a/b+c/d = (a*d+b*c)/(b*d) = (a*(d/g)+b*(c/g))/(b*(d/g)) using gcd g
    %G=gcd(D1,D2);
    G=mygcd(D1,D2);
    D2G=D2/G;
    K=K1+K2;
    N=N1*D2G+N2*(D1/G);
    D=D1*D2G;
    FRAC=fr(K,N,D); % let fr handle reduction
    
  elseif isfinite(FRAC1) || isnan(FRAC2)
    % adding real number to inf or any to nan doesn't change result
    FRAC=FRAC2;
    
  elseif isfinite(FRAC2) || isnan(FRAC1)
    % adding real number to inf or any to nan doesn't change result
    FRAC=FRAC1;
    
  elseif sign(FRAC1)==sign(FRAC2),
    % adding two infinities of same sign:
    % inf+inf = inf or (-inf)+(-inf) = -inf
    FRAC=FRAC1;
  else
    % inf+(-inf) = nan = (-inf)+fin
    FRAC=fr(0,0,0); % 0/0
    
  end; % scalar case
  
elseif (n1 == 1) && (n2 > 1)
  % FRAC1 was a scalar but not FRAC2.
  % Scalar expansion on FRAC1
  FRAC = FRAC2;
  for i = 1:n2
    FRAC(i) = FRAC1 + FRAC2(i);
  end;
  
elseif (n1 > 1) && (n2 == 1)
  % FRAC2 was a scalar but not FRAC1.
  % Scalar expansion on FRAC2
  FRAC = FRAC1;
  for i = 1:n1
    FRAC(i) = FRAC1(i) + FRAC2;
  end;
  
elseif isempty(FRAC1) || isempty(FRAC2)
  % empty propagates
  FRAC = fr([]);
  return
else
  % two non-scalar, non-empty arrays. Are they
  % compatible in size for addition?
  s1 = size(FRAC1);
  s2 = size(FRAC2);
  
  if ~isequal(s1,s2)
    error('FRAC1 and FRAC2 are not compatible in size for addition')
  end
  
  FRAC = FRAC1;
  for i = 1:n2
    FRAC(i) = FRAC1(i) + FRAC2(i);
  end
end

end % main function plus(...)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function g=mygcd(a,b)
% simple gcd function avoids overhead due to arg checks etc
a=abs(a);
b=abs(b);
while a~=0 && b~=0
  if a>=b
    a=mod(a,b);
  else
    b=mod(b,a);
  end;
end;
if a==0
  g=b;
else
  g=a;
end;

end % helper function mygcd(...)
