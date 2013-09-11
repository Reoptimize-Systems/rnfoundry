function [dig,rem] = digits(FRAC,n,d)
% fr/digits: extract the first n digits of the remainder part of a fraction
% usage: dig = digits(FRAC)
%        dig = digits(FRAC,n)
%        dig = digits(FRAC,n,d)
%        [dig,rem] = digits(...)
%        [dig,rep] = digits(...)
%
% arguments: (input)
%  FRAC - a scalar fraction object
%
%  n - a positive integer or empty array, specifying the number of digits
%      to determine.  If empty or not provided, returns the digits up to
%      the first repeat.
%
%  d - a positive integer or fraction-compatible object, specifying the
%      digit base to use (default 10)
%
% arguments: (output)
%  dig - an array of length n containing the sequence of digits of
%        abs(FRAC-fix(FRAC)) up to but not including the first repeat
%  
%  rem/rep - if n is specified, returns the remainder as a fraction,
%            otherwise returns the repeating sequence of digits
%
%
% Example:
%  f = fr(1,3)
%  [dig,rep]=digits(f)
%
%  dig =
%     3
%
%  rep =
%     3
%
%
%  [dig,rem]=digits(f,5)
% 
%  dig =
% 
%      3     3     3     3     3
% 
%  rem =
%    1 / 300000
%
%
%  [dig,rem]=digits(fr(1,3),5,3) % base 3 expansion of 1/3 is 0.1
%
%  dig =
%
%      1     0     0     0     0
%
%  rem =
%     0
%
%
%  See also: fr

% Author: Ben Petschel 27/7/09
%
% Version history:
%   27/7/09 - first release


if (nargin < 1) || (nargin > 3)
  error('fr:digits:nargin','digits requires between 1 and 3 arguments')
end

if nargin == 1

  findrepeat = true;
  n = [];
  d = 10;
  
elseif nargin == 2

  findrepeat = isempty(n);
  d = 10;
  
else
  
  findrepeat = isempty(n);
  
end;

if isempty(FRAC)
  dig=[];
  rem=[];
  return;
end;
if ~isa(FRAC,'fr')
  error('fr:digits:inputclass','input to digits must be a fraction');
end;
if numel(FRAC)>1
  error('fr:digits:inputsize','input to digits must be a scalar fraction');
end;
if ~findrepeat && (numel(n) > 1 || n<0 || floor(n)~=n)
  error('fr:digits:inputsize','n must be a positive integer or an empty array');
end;
if numel(d) ~= 1
  error('fr:digits:inputsize','d must be a scalar object');
end;

try
  dposint = (d>0) && (d==floor(d));
catch
  dposint = false;
end;
if ~dposint
  % let off with a warning (curious to see how non-integer bases are treated)
  warning('fr:digits:inputtype','d is not a positive integer')
end;

if ~isfinite(FRAC)
  % don't determine digits for infinities or nan's
  dig = nan;
  if findrepeat
    rem = nan;
  else
    rem = FRAC;
  end;
  return;
end;

% now find the remainder and multiply by d until done
rem = abs(FRAC-fix(FRAC));
k=0; % number of digits found so far

dig=FRAC.whole+FRAC.numer+FRAC.denom; % this will be overwritten later; just makes dig the same common type as its parts

if ~findrepeat
  dig=repmat(dig,1,n);
else
  % dig(1) will be overwritten later
  dig(1)=[];
  remlist = fr([]);
end;

keepgoing = findrepeat||k<n;
while keepgoing
  if findrepeat
    remlist=[remlist,rem]; % add remainder to list of those found so far
  end;
  k=k+1;
  rem = d*rem;
  dig(k) = rem.whole;
  rem = rem-fix(rem); % remove the whole part
  if findrepeat
    keepgoing = (rem~=0)&&~any(remlist==rem); % see if remainder has already occurred
  else
    keepgoing = k<n;
  end;
end;

if findrepeat
  % return the list of repeating digits
  if rem == 0
    rem = 0*dig(1); % repeating digit is zero (use same class as dig(1))
  else
    ind = find(remlist==rem);
    rem = dig(ind:end);
  end;
else
  % return the remainder
  rem = rem / (d^n);
end;


end
