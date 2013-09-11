function FRAC = mtimes(FRAC1,FRAC2)
% fr/mtimes: Matrix multiplication of two fraction objects
% usage: FRAC = FRAC1*FRAC2;
% usage: FRAC = mtimes(FRAC1,FRAC2);
% 
% arguments:
%  FRAC,FRAC1,FRAC2 - fraction objects, arrays, or numeric variables,
%     FRAC1 and FRAC2 can be arrays if the sizes are valid for matrix mult
%
%
%  See also: times

% Author: Ben Petschel 25/7/09
%
% Version history:
%   25/7/09 - first release (using vpi/mtimes as a template)


if nargin ~= 2
  error('2 arguments required. mtimes is a dyadic operator.')
end;

% make sure both FRAC1 and FRAC2 are fraction objects
if ~isa(FRAC1,'fr')
  FRAC1 = fr(FRAC1);
end
if ~isa(FRAC2,'fr')
  FRAC2 = fr(FRAC2);
end

n1 = numel(FRAC1);
n2 = numel(FRAC2);
if (n1 == 1) || (n2 == 1)
  % scalar inputs or scalar expansion,
  % implement as just a call to times
  FRAC = FRAC1.*FRAC2;
elseif (n1*n2) == 0
  % empty will propagate
  FRAC = fr([]);
else
  % an array multiplication
  s1 = size(FRAC1);
  s2 = size(FRAC2);
  
  % do not allow multiplication of arrays
  % in higher dimensions to be consistent
  % with standard matlab syntax.
  if (length(s1) > 2) || (length(s2) > 2)
    error('fr:mtimes:multidim','Matrix multiplication of multidimensional arrays is not supported');
  end
  
  % Do the arrays conform in size?
  if s1(2) ~= s2(1)
    error('fr:mtimes:sizes','The arrays do not conform in size for matrix multiplication')
  end
  
  % preallocate the result as zeros
  FRAC = repmat(fr(0),s1(1),s2(2));
  % just loops now
  for i = 1:s1(1)
    for j = 1:s2(2)
      for k = 1:s1(2)
        FRAC(i,j) = FRAC(i,j) + FRAC1(i,k).*FRAC2(k,j);
      end
    end
  end
end


