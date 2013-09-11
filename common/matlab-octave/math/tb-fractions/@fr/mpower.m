function b = mpower(a,k)
% fr/power: raises a fraction object to a positive integer power
% usage: b = a^k;
% usage: b = mpower(a,k);
% 
% arguments: (input)
%  a - a fraction object scalar or square array
%
%  k - positive integer scalar,
%      k may be no larger than 2^53 - 1
%
% arguments: (output)
%  b - a fraction object, representing the power operation a^k
%
%  See also: power

% Author: Ben Petschel 25/7/09
%
% Version history:
%   25/7/09 - first release (using vpi/mpower as a template)


na = numel(a);
nk = numel(k);

if (na == 0) || (nk == 0)
  % propagate empty
  b = fr([]);
elseif (na == 1) && (nk == 1)
  % scalar ops, use power
  b = a.^k;
elseif (na > 1) && (nk == 1) && (k >= 0) && ...
    (ndims(a) == 2) && (k <= (2^53 - 1))
  % we can raise a square matrix to a positive
  % integer power
  if (diff(size(a)) ~= 0)
    error('fr:mpower:inputsize','Cannot raise a non-square matrix to an integer power')
  end
  
  if isnumeric(k)
    kbin = dec2bin(k) - '0';
  elseif isa(k,'vpi')
    kbin = vpi2bin(k) - '0';
  else
    error('fr:mpower:ktype','k is not a valid type for calculating a^k');
  end
  kbin = fliplr(kbin);
  
  % what size square matrix is a?
  n = size(a,1);
  if kbin(1)
    b = a;
  else
    b = fr(eye(n));
  end
  a2 = a;
  for i = 2:length(kbin)
    % repeatedly square a
    a2 = a2*a2;
    
    % do we multiply this term in?
    if kbin(i)
      b = b*a2;
    end
  end
else
  % not supported
  error('Not supported: raising a scalar or matrix to a non-scalar power')
end


