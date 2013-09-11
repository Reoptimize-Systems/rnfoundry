function Asum = sum(A,dim)
% fr/sum: sum of a fraction object array
% usage: Asum = sum(A,dim);
% 
% Arguments:
%  A - a fraction object array
%
%  dim - (optional) dimension of A to sum over
%      DEFAULT: dim is the first non-singleton
%      dimension of A.
%
% Arguments: (output)
%  Asum - the sum fraction object array
% 
% Example:
%  A = fr([1:3],4); % 1/4, 2/4, 3/4
%  sum(A)
%
%  ans =
%     1 + 1 / 2
%
%  See also: prod, cumsum

% Author: Ben Petschel 27/7/09
%
% Version history:
%   27/7/09 - first release (using vpi/sum as a template)

if (nargin<1) || (nargin>2)
  error('fr:sum:nargin','sum takes one or two arguments')
end

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
  indsum=ind;
  indsum{dim}=1;
  Asum=A(indsum{:}); % start off with A(:,...,:,1,:,...,:)
  for i=2:sA(dim)
    ind{dim}=i;
    Asum(indsum{:})=Asum(indsum{:})+A(ind{:});
  end;
  
end;

end
