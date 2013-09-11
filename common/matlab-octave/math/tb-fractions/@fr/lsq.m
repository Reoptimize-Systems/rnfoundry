function [X,N]=lsq(A,B)
% fr/lsq - finds the least-squares solution to A*X=B where A and B are fraction arrays
% usage: X=lsq(A,B)
%        [X,N]=lsq(A,B)
%
% lsq(A,B) solves the system (A'*A)*X=(A'*B)
%
% Example: system x=0, x=1
%   lsq(fr([1;1]),[0;1]) % returns 1/2
%
% see also mldivide

% Author: Ben Petschel 17/8/2009
% Version history:
%  17/8/2009 - first release


[na,ma]=size(A);
[nb,mb]=size(B);

if na~=nb
  error('fr/lsq:inputsize','number of input rows does not match');
end;
if ~isa(A,'fr')
  A=fr(A);
end;
if ~isa(B,'fr');
  B=fr(B);
end;

[X,N]=(A'*A)\(A'*B);

end

