function Aprod = prod(A,dim)
% fr/prod: prod of a fraction array
% usage: Aprod = prod(A,dim);
% 
% Arguments:
%  A - a fraction object array
%
%  dim - (optional) dimension of A to sum over
%      DEFAULT: dim is the first non-singleton
%      dimension of A.
%
% Arguments: (output)
%  Aprod - the product fraction object array
% 
% Example:
%  A = fr([1,2,3],[4,5,6]); % 1/4, 2/5, 3/6
%  prod(A)
%
%  ans =
%     1 / 20
%
%  See also: length, numel, size, repmat, sum

% Author: Ben Petschel 27/7/09
%
% Version history:
%   27/7/09 - first release (using vpi/prod as a template)


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
  indp{dim}=1;
  Aprod=A(indp{:}); % start off with A(:,...,:,1,:,...,:)
  for i=2:sA(dim)
    ind{dim}=i;
    Aprod(indp{:})=Aprod(indp{:}).*A(ind{:});
  end;
  
end;

end
