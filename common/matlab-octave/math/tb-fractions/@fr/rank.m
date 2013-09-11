function r=rank(A)
% fr/rank - determines rank of a fraction array
% usage: r=rank(A)
%
% r=rank(A) returns the number of linearly independent rows in the fraction
% array A.  The results are exact up to any limitations of the underlying
% data types of the fraction numerator/denominators.  For larger arrays and
% arrays with large denominators, intermediate calculations may exceed
% these limitations and use of vpi fractions is recommended.
%
% Returns NaN if any elements of rref(A) are not finite.
%
% see also fr, vpi

% Author: Ben Petschel 15/12/2009
% Version history:
%  15/12/2009 - first release


% simple algorithm here uses rref (ok for medium-size problems)
% counts the number of leading non-zeros in the rows of the rref
B=rref(A);
if any(~isfinite(B))
  r=nan;
else
  [mb,nb]=size(B);
  r=0;
  i=1;
  j=1;
  while i<=mb && j<=nb
    if B(i,j)==0
      j=j+1;
    else
      r=r+1;
      i=i+1;
      j=j+1;
    end;
  end; % while
end; % if ... else ...
