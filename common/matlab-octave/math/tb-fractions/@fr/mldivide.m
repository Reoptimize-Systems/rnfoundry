function [X,N]=mldivide(A,B)
% fr/mldivide - solves A*X=B where A and B are fraction arrays
% usage: X=mldivide(A,B)
%        X=A\B
%        [X,N]=mldivide(A,B)
%        [X,N]=A\B
%
% A\B, where A and B are fraction arrays with the same number of rows,
% solves the linear system A*X=B by gaussian elimination with row and
% column pivoting.
%
% Treatment of singular and non-square arrays differs from that of the
% built-in "\".  Non-square and singular systems are solved by looking for
% a particular solution and returning a warning if the solution does not
% exist or is not unique.  To find a least-squares solution, use LSQ.
%
% If A is non-square or singular and the system A*X=B is inconsistent, then
% there is no solution and the empty array is returned.
%
% If A is non-square or singular and the system A*X=B is indeterminate:
% (i)  [X,N] = A\B returns a particular solution X and the null space N
%      such that X+N*V is a solution where V is any fraction array of
%      compatible size.
% (ii) X = A\B returns NaN values in some of the elements of X to indicate
%      that the solution was not unique.
%
% A non-square matrix may still yield a unique solution, e.g. if A and B
% are padded with additional rows of zeros.
%
% If A is scalar, A\B returns B/A (scalar division).
%
% To switch off the warnings that can occur when A is singular:
%   warning('off','fr:mldivide:singular')
%
%
% Example (system x+y=1, 0*x+0*y=0 with solution x=1-y):
%   [1,1]\1             % built-in "\" returns [1;0], no warning
%   [1,1;0,0]\[1;0]     % built-in "\" returns [NaN;NaN], with warning
%   fr([1,1])\1         % returns [1;NaN], with warning
%   fr([1,1;0,0])\[1;0] % returns same result
%   [X,N]=fr([1,1;0,0])\[1;0] % returns X=[1;0] and N=[-1;1], no warning
%
%
% The results are exact up to any limitations of the underlying data types
% of the fraction numerator/denominators.  For larger arrays and arrays
% with large denominators, intermediate calculations may exceed these
% limitations and use of vpi fractions is recommended.
%
% see also fr, vpi

% Author: Ben Petschel 17/8/2009
% Version history:
%  17/8/2009 - first release

if nargout==2
  findparticular=true; % indicates that a particular solution should be found
else
  if numel(A)==1
    X=B/A; % scalar division
    N=[];
    return;
  end;
  findparticular=false; % returns nan values if more than one solution exists
end;

[na,ma]=size(A);
[nb,mb]=size(B);

if na~=nb
  error('fr/mldivide:inputsize','number of input rows does not match');
end;
if ~isa(A,'fr')
  A=fr(A);
end;
if ~isa(B,'fr');
  B=fr(B);
end;


R=B;
U=A;


% now do gaussian elimination with pivoting to get  U*(P^(-1))*x = R
i=1;
p=1:ma; % indexes of permutation such that P*v=v(p)
pind=1:ma; % pind(i) is position of i in p.
nosol=false; % set to true if any inconsistent equations found
imax=min(na,ma);
while i<=imax
  % see if pivoting is required
  Uiinz=true;
  if U(i,i)==0
    % look for a row pivot
    ind=i+find(U(i+1:end,i)~=0,1,'first'); % returns [] if no such value
    if ~isempty(ind)
      % pivot by swapping rows i and ind in U and R
      tmp=U(i,:);
      U(i,:)=U(ind,:);
      U(ind,:)=tmp;
      tmp=R(i,:);
      R(i,:)=R(ind,:);
      R(ind,:)=tmp;
    else
      % look for a column pivot
      ind=i+find(U(i,i+1:end)~=0,1,'first'); % returns [] if no such value
      if isempty(ind)
        Uiinz=false;
      else
        % pivot by swapping columns i and ind; then swap positions of i and
        % ind in p
        tmp=U(:,i);
        U(:,i)=U(:,ind);
        U(:,ind)=tmp;
        p1=pind(i);
        p2=pind(ind);
        p(p1)=ind;
        p(p2)=i;
        pind(i)=p2;
        pind(ind)=p1;
      end;
    end;
  end;
  if ~Uiinz && any(R(i,:)~=0)
    % have inconsistent equations
    nosol=true;
  elseif Uiinz
    % row reduce on U(i,i)
    Uiend=U(i,i+1:end)/U(i,i);
    Ri=R(i,:)/U(i,i);
    for j=1:na
      if j~=i
        Uji=U(j,i);
        if U(j,i)~=0
          % only work on elements i:end of row j (all others are zero)
          U(j,i+1:end)=U(j,i+1:end)-Uji*Uiend;
          U(j,i)=fr(0);
        end;
        % reduce R
        R(j,:)=R(j,:)-Uji*Ri;
      end;
    end;
    U(i,i+1:end)=Uiend;
    U(i,i)=1;
    R(i,:)=Ri;
  end;
  i=i+1;
end;

% have U*(P^(-1)*x)=R where U is upper triangular with all zeros or ones
% on diagonal, so solve y = U\R + N*z
N=repmat(fr([]),ma,1); % ma x 0 fraction array
y=R;
hasinf=false;
hasnan=false;

% check that all spare rows of R (if any) are zero
if any(any(R(ma+1:end,:)~=0))
  hasinf=true;
  nosol=true;
end;

j=ma;
while j>0
  if j>na || U(j,j)==0
    % add to null space
    % find places where U(i,j)~=0; have X(i)+U(i,j)*X(j)=0
    Nj=fr(zeros(ma,1));
    ind=find(U(:,j)~=0);
    Nj(ind)=-U(ind,j);
    Nj(j)=1;
    N=[N,Nj];
    if j>na
      hasnan=true;
      if findparticular
        y(j,:)=0;
      else
        % put nan values in solution
        y(j,:)=fr(nan);
      end
    else
      if findparticular
        % leave y(j,:) as zero, to get particular solution
        yjnz=y(j,:)~=0; % nonzero positions
        hasinf=hasinf | any(yjnz);
        y(j,yjnz)=fr(inf);
      else
        yjnz=y(j,:)~=0;
        hasinf=hasinf | any(yjnz);
        hasnan=hasnan | any(~yjnz);
        y(j,:)=abs(y(j,:))/0; % get nan or inf
      end;
    end;
  else
    % U(i,i)=1 and U(j,i)=0 for j~=i so leave y alone
    %yi=y(i,:)/U(i,i);
    %y(i,:)=yi;
    %for j=1:i-1
    %  y(j,:)=y(j,:)-U(j,i)*yi;
    %end;
  end;
  j=j-1;
end;

if hasinf && hasnan
  warning('fr:mldivide:singular','Coefficient matrix is singular and the system has both indeterminate and inconsistent components. Recommend using [X,N]=lsq(A,B)');
elseif hasinf
  warning('fr:mldivide:singular','Coefficient matrix is singular and the system has inconsistent components.  Recommend using least squares X=lsq(A,B)');
elseif ~findparticular&&hasnan
  % only issue warning if null space was not computed
  warning('fr:mldivide:singular','Coefficient matrix is singular and the system has indeterminate components.  Recommend using [X,N]=A\\B');
end;

if nosol
  X=zeros(ma,0); % ma x 0 empty matrix
else
  X=y(p,:); % permute rows of X
end;
if ~isempty(N)
  N=N(p,:); % permute rows of N
end;

end

