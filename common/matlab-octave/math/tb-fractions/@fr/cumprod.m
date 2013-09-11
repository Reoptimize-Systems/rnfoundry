function Aprod = cumprod(A,dim)
% fr/cumprod: cumulative product of a fraction array
% usage: Aprod = cumprod(A,dim);
% 
% Arguments:
%  A - a fraction object array
%
%  dim - (optional) dimension of A to prod over
%      DEFAULT: dim is the first non-singleton
%      dimension of A.
%
% Arguments: (output)
%  Aprod - the product fraction object array
% 
%
%  See also: cumsum, prod, sum

% Author: Ben Petschel 27/7/09
%
% Version history:
%   27/7/09 - first release (using vpi/cumprod as a template)


if (nargin<1) || (nargin>2)
  error('prod takes one or two arguments')
end;

if numel(A) == 1
  % the product of a scalar is a no-op
  Aprod = A;
else
  % a vector or array
  
  % default for dim?
  sA=size(A);
  if (nargin==1) || isempty(dim)
    dim = find(sA>1,1,'first');
    if isempty(dim)
      dim = 1;
    end;
  end;
  
  % prod over the dimension dim, using string-based array referencing
  ind=repmat({':'},1,numel(sA));
  indp=ind;
  Aprod=A;
  for i=2:sA(dim)
    indp{dim}=i-1; % previous index
    ind{dim}=i; % current index
    Aprod(ind{:})=Aprod(indp{:}).*A(ind{:});
  end;
  
end;

end
