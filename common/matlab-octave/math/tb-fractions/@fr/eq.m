function B = eq(FRAC1,FRAC2)
% fr/eq: test for equality between a pair of fraction objects
% usage: B = (FRAC1 == FRAC2)
% usage: B = eq(FRAC1,FRAC2)
% 
% arguments: (input)
%  FRAC1, FRAC2 - fraction objects or scalar numeric values,
%                 can be arrays provided the sizes are compatible
%
% arguments: (output)
%  B  - logical variable, true when the two inputs represent the same value
%       (treating pairs of NaNs as equal).
% 
% Example:
%  FRAC = fr(1,2);
%  FRAC == 1/2
%  ans =
%     1
%
%  FRAC == FRAC
%  ans =
%     1
%
%
%  See also: ge, gt, lt, le

% Author: Ben Petschel 25/7/09
%
% Version history:
%   25/7/09 - first release (using vpi/eq as a template)



if nargin~=2
  error('== is a dyadic operator, exactly 2 arguments are required')
end

% convert objects to fractions
if ~isa(FRAC1,'fr')
  FRAC1=fr(FRAC1);
end;
if ~isa(FRAC2,'fr')
  FRAC2=fr(FRAC2);
end;

n1=numel(FRAC1);
n2=numel(FRAC2);

if n1==1 && n2==1
  % scalar vs scalar
  s=[1,1];
elseif n1==1
  % array vs scalar
  s=size(FRAC2);
elseif n2==1
  % scalar vs array
  s=size(FRAC1);
else
  s=size(FRAC1);
  if ~isequal(s,size(FRAC2)),
    error('fr:eq:inputsize','array input sizes must match');
  end;
end;

% treat nan's as equal

if n1==1 && n2==1
  % use shortcut for scalar case
  B = (FRAC1.whole==FRAC2.whole) && (FRAC1.numer==FRAC2.numer) ...
    && (FRAC1.denom==FRAC2.denom);
  
else
  % array case
  B = false(s);
  B(:) = ([FRAC1(:).whole]==[FRAC2(:).whole]) ...
    & ([FRAC1(:).numer]==[FRAC2(:).numer]) ...
    & ([FRAC1(:).denom]==[FRAC2(:).denom]);
end;

end
