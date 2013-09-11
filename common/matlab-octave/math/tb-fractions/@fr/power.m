function FRAC = power(FRAC1,k)
% fr/power: raises a fraction object to a positive integer (numeric) power
% usage: FRAC = FRAC1.^k;
% usage: FRAC = power(FRAC1,k);
%
% arguments: (input)
%  FRAC1 - a fraction object
%
%  k   - positive integer scalar, k must be no larger than 2^53
%
% arguments: (output)
%  FRAC  - a fraction object
%
%
%  See also: mpower

% Author: Ben Petschel 25/7/09
%
% Version history:
%   25/7/09 - first release (using vpi/power as a template)
%   20/8/09 - general speedup by removing special cases


% is k an integer?
if any(k<0 | floor(k)~=k | k>=2^53-1)
  error('fr:power:nonneg','k must be a nonnegative integer between 0 and 2^53-1');
end;

% if ~isa(FRAC1,'fr')
%   % Make sure FRAC1 is a fraction object.
%   FRAC1 = fr(FRAC1);
% end;

nf=numel(FRAC1);
nk=numel(k);

% There are several possibilities to worry about
if (nf == 1) && (nk == 1)
  % a pair of scalars to work with

% % removed special cases because of the overhead for comparing to 0,1,-1
%   % special cases, is FRAC1 0, 1 or -1?
%   if FRAC1 == 0
%     % 0^k == 0, except 0^0 = NaN
%     if k ~= 0
%       FRAC = fr(0);
%     else
%       FRAC = fr(NaN);
%     end
%     return
%   elseif FRAC1 == 1
%     % 1^k == 1
%     FRAC = fr(1);
%     return
%   elseif FRAC1 == -1
%     % (-1)^k
%     if rem(k,2) == 0
%       FRAC = fr(1);
%     else
%       FRAC = fr(-1);
%     end
%     return
%   end

  % other special cases - where k is zero or one
  if k == 0
    % anything to the zero power is 1
    FRAC = fr(1,0,1);
    return
 
  elseif k == 1
    % anything to the 1 power is itself.
    FRAC = FRAC1;
    return

  elseif k == 2
    % simple powers are easy to do as a multiply
    FRAC = FRAC1.*FRAC1;
    return

  elseif k == 3
    % simple powers are easy to do as a multiply
    FRAC = FRAC1.*FRAC1.*FRAC1;
    return

  elseif k == 4
    % simple powers are easy to do as a multiply
    FRAC = FRAC1*FRAC1;
    FRAC = FRAC*FRAC;
    return

  elseif k == 8
    % simple powers are easy to do as a multiply
    FRAC = FRAC1*FRAC1;
    FRAC = FRAC*FRAC;
    FRAC = FRAC*FRAC;
    return

  end

  % otherwise (K+N/D)^k = (K*D+N)^k / D^k
  K=FRAC1.whole;
  N=FRAC1.numer;
  D=FRAC1.denom;
  FRAC=fr(0,(K*D+N).^k,D.^k);

elseif nf==0 || nk==0
  % empty propagates to empty
  FRAC = fr([]);
  
elseif (nf == 1) && (nk > 1)
  % A scalar raised to an array power as
  % an elementwise operation. Do scalar
  % expansion for FRAC1.
  FRAC = repmat(fr(0),size(k));
  for i = 1:nk
    FRAC(i) = FRAC1.^k(i);
  end;

elseif (nf > 1) && (nk == 1)
  % A fraction array raised to a scalar power as
  % an elementwise operation. Do scalar
  % expansion for k.
  FRAC = FRAC1;
  for i = 1:nf
    FRAC(i) = FRAC1(i).^k;
  end;

elseif (nf > 1) && (nk > 1)
  % A fraction array raised to an array power as
  % an elementwise operation.

  % first verify the arrays conform in size
  sf = size(FRAC1);
  sk = size(k);
  if ~isequal(sf,sk)
    error('fr:power:argsize','FRAC1 and k do not conform in size for elementwise power op');
  end;

  FRAC = FRAC1;
  for i = 1:nk
    FRAC(i) = FRAC1(i).^k(i);
  end;

end;

end




