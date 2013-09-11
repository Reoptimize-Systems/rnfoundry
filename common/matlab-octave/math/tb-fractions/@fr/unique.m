function [S,ind]=unique(A,varargin)
% fr/UNIQUE - find the unique elements of a fraction array
%
% Usage: S=unique(A)
%        S=unique(A,'rows')
%        S=unique(...,'first')
%        S=unique(...,'last')
%        [S,ind]=unique(...)
%        
% Intended for use with fraction object arrays however it should work with
% any object types with "sort", "sortrows", "isequal" and "==" defined.
%
% See the help for the built-in unique function.

% Author:
%   Ben Petschel 20/7/2009
% Change history:
%   20/7/2009: initial release (donated to VPI toolbox)
%   28/7/2009: updated documentation for fraction object types

byrows=false;
firstind=true;
if nargin>=2,
  for i=1:nargin-1,
    switch varargin{i}
      case 'rows'
        byrows=true;
      case 'first'
        firstind=true;
      case 'last'
        firstind=false;
    end;
  end;
end;

if isempty(A),
  S=A;
  ind=[];
  return;
end;

if byrows,
  % remove identical rows
  [B,ind]=sortrows(A);
  i=1;
  for j=2:size(B,1),
    if isequal(B(i,:),B(j,:)),
      if firstind,
        ind(j)=0;
      else
        ind(i)=0;
        i=j;
      end;
    else
      i=j;
    end;
  end;
  ind=ind(ind>0);
  S=A(ind,:);
else
  if ~isvector(A),
    [B,ind]=sort(A(:));
  else
    [B,ind]=sort(A);
  end;
  i=1;
  for j=2:numel(B),
    if B(i)==B(j),
      if firstind,
        ind(j)=0;
      else
        ind(i)=0;
        i=j;
      end;
    else
      i=j;
    end;
  end;
  ind=ind(ind>0);
  S=A(ind);
end;


end % function unique(...)

