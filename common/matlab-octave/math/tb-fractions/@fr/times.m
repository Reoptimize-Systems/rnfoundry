function FRAC = times(FRAC1,FRAC2)
% fr/times: Multiplies two fraction objects
% usage: FRAC = FRAC1.*FRAC2;
% usage: FRAC = times(FRAC1,FRAC2);
% 
% arguments: (input)
%  FRAC1,FRAC2 - fraction objects, or any compatible objects,
%                can be arrays if the sizes are compatible
%
% arguments: (output)
%  FRAC - a fraction that represents the product FRAC1.*FRAC2
%
%
%  See also: mtimes

% Author: Ben Petschel 25/7/09
%
% Version history:
%   25/7/09 - first release (using vpi/times as a template)
%   20/8/09 - improved speed by rearranging multiplication

if nargin ~= 2
  error('2 arguments required. times is a dyadic operator.')
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
  % a pair of scalar elements to multiply

  % check for inf and nan
  % nan*x=nan, inf*0=nan, inf*x=sign(x)*inf
  if isfinite(FRAC1) && isfinite(FRAC2)
    % (k1+n1/d1)*(k2+n2/d2)
    % = k1*k2 + (k1*n2)/d2 + (k2*n1)/d1 + (n1*n2)/(d1*d2)
    % let k1*n2 = c1*d2+r1
    % let k2*n1 = c2*d1+r1
    % -> k1*k2 + c1 + c2 + r1/d2 + r2/d1 + (n1*n2)/(d1*d2)
    % -> k1*k2 + c1 + c2 + (r1*d1 + r2*d2 + n1*n2)/(d1*d2)
    K1=FRAC1.whole;
    K2=FRAC2.whole;
    N1=FRAC1.numer;
    N2=FRAC2.numer;
    D1=FRAC1.denom;
    D2=FRAC2.denom;
    k1n2=K1*N2;
    k2n1=K2*N1;
    r1=mod(k1n2,D2);
    c1=(k1n2-r1)/D2;
    r2=mod(k2n1,D1);
    c2=(k2n1-r2)/D1;
    FRAC=fr(K1*K2+c1+c2,r1*D1+r2*D2+N1*N2,D1*D2);
    % %old method: use gcd to avoid overflow if possible, then pass to plus
    %gk1d2=gcd(K1,D2);
    %gk2d1=gcd(K2,D1);
    %gn1d2=gcd(N1,D2);
    %gn2d1=gcd(N2,D1);
    %FRAC=fr(0,K1*K2,1) ...
    %  + fr(0,(K1/gk1d2)*N2,(D2/gk1d2)) ...
    %  + fr(0,(K2/gk2d1)*N1,(D1/gk2d1)) ...
    %  + fr(0,(N2/gn2d1)*(N1/gn1d2),(D1/gn2d1)*(D2/gn1d2));
    
  elseif isnan(FRAC1)
    % nan*x=nan
    FRAC=FRAC1;
  elseif isnan(FRAC2)
    % x*nan=nan
    FRAC=FRAC2;

    % otherwise left with one or two inf values
  elseif iszero(FRAC1) || iszero(FRAC2)
    % inf * 0 = nan
    FRAC=fr(0,0,0);
  elseif isfinite(FRAC1)
    % x * inf = sign(x)*inf
    if sign(FRAC1)>0,
      FRAC=FRAC2;
    else
      FRAC=-FRAC2;
    end;
  else
    % inf*x = sign(x)*inf
    if sign(FRAC2)>0,
      FRAC=FRAC1;
    else
      FRAC=-FRAC1;
    end;
  end; % scalar case
    
elseif (n1 == 0) || (n2 == 0)
  % empty propagates
  FRAC = fr([]);
  return
elseif (n1 == 1) && (n2 > 1)
  % scalar expansion for FRAC1
  FRAC = fr(FRAC2);
  for i = 1:n2
    FRAC(i) = FRAC1.*FRAC2(i);
  end
elseif (n1 > 1) && (n2 == 1) 
  % scalar expansion for FRAC2
  FRAC = fr(FRAC1);
  for i = 1:n1
    FRAC(i) = FRAC1(i).*FRAC2;
  end
else
  % must be two arrays
  
  % do they conform for multiplication?
  if ~isequal(size(FRAC1),size(FRAC2))
    error('The two arrays do not conform in size for elementwise multiplication')
  end
  
  % do the scalar multiplies
  FRAC = fr(FRAC1);
  for i = 1:n1
    FRAC(i) = FRAC1(i).*FRAC2(i);
  end
end




