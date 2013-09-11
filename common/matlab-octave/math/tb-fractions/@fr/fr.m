function FRAC = fr(varargin)
% fr: Creates a fraction object
% usage: F = fr
%        F = fr(K)     % K is a numeric variable
%        F = fr(N,D)   % N and D are numeric variables
%        F = fr(K,N,D) % K, N and D are numeric variables
%        F = fr(...)   % K, N and D are compatible objects (see below)
%
%
% Arguments:
%  K - the whole number part of the fraction
%  N - the numerator of the fraction
%  D - the denominator of the fraction
%
%  fr(K,N,D) returns a fraction object representing K + N/D.
%  K, N, D can be integers (double, single, int*) or any "compatible"
%  objects for which certain arithmetic and comparison operations are
%  defined (see below).  K, N, D are assumed to be 0, 0, 1 if not provided.
%
%
%  The fraction will be automatically reduced to lowest common form, N>0.
%  - non-integer floats N (single/double) are converted to fractions using
%    RAT however this has a default tolerance of 1e-6*abs(N).
%    For better accuracy it is recommended that N and D are precomputed
%    with rats using a smaller tolerance.
%  - Integers N are represented by N + 0 / 1
%  - Negative fractions -N/D are represented by -1 + (D-N) / D (if N<D)
%  - [-inf,nan,inf] are represented by 0 + [-1,0,1] / 0
%
%  The following object types are known to be compatible:
%    vpi (John D'Errico's Variable Precision Integer toolbox,
%         available on MATLAB File Exchange)
%
%  Non-standard objects must include 0, 1, -1 and require the following
%  operations to be defined in order to create a fraction object:
%    gcd
%    rem
%    sign
%    abs
%    +, - , .*, ./
%    ==, <, <=, >, >=, ~=
%  
%  The following additional operation definitions are recommended:
%    *, .^
%    sort
%    floor
%    factor
%    gcd (3-output form)
%    rat (if floor(x) or mod(x,1) is not always equal to x)
%
%  Examples:
%  f = fr      % creates a zero fraction (0+0/1)
%  f = fr([])  % creates an empty fraction object
%  f = fr(1,3) % creates a fraction representing 1/3 (0+1/3)
%  f = fr(1,vpi(3)^100) % creates a fraction representing 1/3^100
%                       % (requires VPI toolbox)
%
%  See also: double, single

% Author: Ben Petschel 25/7/2009
%
% Version history:
%   25/7/2009 - first release (using vpi/vpi as a template)
%   28/7/2009 - added checks to ensure fraction parts are all the same class
%   6/10/2012 - substantial speedup by using "superiorto" only on first call

% process the input arguments
if nargin==0
  % creates a zero variable
  K = 0;
  N = 0;
  D = 1;
elseif nargin==1
  K=varargin{1};
  if isa(K,'fr'),
    % already a fraction object
    FRAC=K;
    return;
  else
    % convert to a fraction later
    N=0;
    D=1;
  end;
elseif nargin==2
  K=0;
  N=varargin{1};
  D=varargin{2};
elseif nargin==3
  K=varargin{1};
  N=varargin{2};
  D=varargin{3};
else
  error('fr:nargin','fr cannot have more than 3 input arguments');
end;


% create the fraction object
isfrac=false;
isreduced=true;
if isempty(N) || isempty(D) || isempty(K)
  % create an empty fraction
  FRAC.numer = 0;
  FRAC.denom = 1;
  FRAC.whole = 0;
  FRAC(1) = [];

elseif any([numel(N),numel(D),numel(K)]>1)
  % create an array: check the sizes agree and then replicate frac objects
  
  % check sizes agree and find common sizes
  sn=size(N);
  sd=size(D);
  sk=size(K);
  scell={sn,sd,sk};
  sf=scell{find(cellfun(@(x)any(x>1),scell),1,'first')}; % get first non-1x1 size
  sok=all(cellfun(@(x)(all(x==1)||isequal(x,sf)),scell)); % check all sizes are 1x1 or sf
  
  if sok,
    % sizes are ok, so proceed
    FRAC.numer = 0;
    FRAC.denom = 1;
    FRAC.whole = 0;
  
    FRAC = repmat(FRAC,sf);
    
    if all(sn==1)
      N=repmat(N,sf);
    end;
    if all(sd==1)
      D=repmat(D,sf);
    end;
    if all(sk==1)
      K=repmat(K,sf);
    end;
    for i=1:prod(sf)
      FRAC(i)=fr(K(i),N(i),D(i));
    end;
    
  else
    error('fr:inputsize','array sizes of N, D and K do not match');
  end;

else
  % All scalars
  
  % FK, FN, FD are only converted to fractions if they are are not integers
  [tfk,FK]=tofrac(K);
  [tfn,FN]=tofrac(N);
  [tfd,FD]=tofrac(D);
  if tfk || tfn || tfd,
    if ~tfn && ~tfd,
      % adding K part was fraction, N and D are integers
      if N==0 && D~=0,
        % N/D = 0
        FRAC = FK;
      else
        FRAC = FK + fr(N,D);
      end;
    else
      % N and/or D are fractions
      FRAC = FK + FN/FD;
    end;
    isfrac=true;
    
  else
    if (isa(N,'double') && (abs(N)>=(2^53))) || (isa(D,'double') && (abs(D)>=(2^53))) || (isa(K,'double') && (abs(K)>=(2^53)))
      warning('fr:fr:largedouble','Any N, D or K that are doubles must be no larger than 2^53 - 1, otherwise roundoff error will result.');
    end
    
    % ensure that N,D,K are the same type
    try
      zeroterm = (N-N)+(D-D)+(K-K);
    catch %#ok<CTCH>
      % types are incompatible
      error('fr:fr:inputclass','Classes of N,D,K are incompatible for addition');
    end;
    
    FRAC.numer = N+zeroterm;
    FRAC.denom = D+zeroterm;
    FRAC.whole = K+zeroterm;

    isreduced=false;
  end;
end

% set the class for this variable
if ~isfrac
  FRAC = class(FRAC,'fr');
end;

% reduce the fraction to lowest form
if ~isreduced
  FRAC = freduce(FRAC);
end;


% make sure that the appropriate fraction method is
% used whenever one of the arguments is a fraction
% and any numeric type as the other argument (need only do this once).
persistent classinitialised % initially empty
if isempty(classinitialised)
  classinitialised = true;
  superiorto('double','single','int8','uint8','int16', ...
    'uint16','int32','uint32','int64','uint64','logical')
  try
    superiorto('vpi')
  catch %#ok<CTCH>
    % in case VPI toolbox has never been used
  end;
end;

end % main function fr(...)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tf,F]=tofrac(X)
% determines whether the input is a whole number or not,
% and convert to a fraction.  tf returns true if it was converted to a
% fraction.

tf=false;

try
  if isnan(X)
    tf=true;
    F=fr(0,0,0);
    return
  elseif isinf(X)
    tf=true;
    F=fr(0,sign(X),0);
    return
  end;
catch %#ok<CTCH>
  % for objects where isnan or isinf is not defined assume value is finite
end;

% otherwise see if X has a remainder mod 1
try
  K=floor(X);
  rem=X-K;
catch %#ok<CTCH>
  % for object types that don't have a "floor" function
  rem=mod(X,1);
  K=X-rem;
end;

if rem==0,
  F=X;
else
  try
    % assume double, otherwise error will result
    [N,D]=rat(rem);
    tf=true;
    F=fr(K,N,D);
  catch %#ok<CTCH>
    warning('fr:convertfraction','input value has a remainder mod 1, however function "rat" is not defined for this object; rounding downwards instead');
    F=X;
  end;
end;

end % helper function tofrac(X)