function Asum = cumsum(A,dim)
% fr/cumsum: cumulative sum of a fraction array
% usage: Asum = cumsum(A,dim)
% 
% Arguments:
%  A - a fraction object array
%
%  dim - (optional) dimension of A to sum over
%      DEFAULT: dim is the first non-singleton
%      dimension of A.
%
% Arguments: (output)
%  Asum - the cumulative sum of the fraction object array
% 
%
%  See also: prod, sum, cumprod

% Author: Ben Petschel 27/7/09
%
% Version history:
%   27/7/09 - first release (using vpi/cumsum as a template)


if (nargin<1) || (nargin>2)
  error('cumsum takes one or two arguments')
end;

if numel(A) == 1
  % the sum of a scalar is a no-op
  Asum = A;
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
  
  % sum over the dimension dim, using string-based array referencing
  ind=repmat({':'},1,numel(sA));
  indp=ind;
  Asum=A; % start off with A(:,...,:,1,:,...,:)
  for i=2:sA(dim)
    indp{dim}=i-1; % previous index
    ind{dim}=i; % current index
    Asum(ind{:})=Asum(indp{:})+A(ind{:});
  end;
  
end










