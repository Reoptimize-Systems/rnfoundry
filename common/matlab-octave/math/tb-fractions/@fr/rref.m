function B=rref(A)
% fr/rref - determines reduced row echelon form of a fraction array
% usage: B=rref(A)
%
% B=rref(A) returns the reduced row echelon form of the fraction array A.
% The results are exact up to any limitations of the underlying data types
% of the fraction numerator/denominators.  For larger arrays and arrays
% with large denominators, intermediate calculations may exceed these
% limitations and use of vpi fractions is recommended.
%
% see also fr, vpi

% Author: Ben Petschel 14/8/2009
% Version history:
%  14/8/2009 - first release

[n,m]=size(A);
B=A;

i=1; % first i-1 rows are in echelon form
j=1; % first j-1 columns are in echelon form

while i<=n && j<=m
  % find a non-zero element to pivot on

  ind=i-1+find(B(i:end,j)~=0,1,'first');
  if ~isempty(ind)
    % row-reduce
    if ind==i
      B(i,:)=B(i,:)/B(i,j);
    else
      % swap row i and ind
      tmp=B(ind,:);
      B(ind,:)=B(i,:);
      % divide row i by B(i,j)
      B(i,:)=tmp/tmp(j);
    end;
    % now reduce all other rows wrt row i
    for ind=[1:i-1,i+1:n]
      if B(ind,j)~=0
        B(ind,:)=B(ind,:)-B(ind,j)*B(i,:);
      end;
    end;
    i=i+1; % go on to next row
  end;
  j=j+1;
end

