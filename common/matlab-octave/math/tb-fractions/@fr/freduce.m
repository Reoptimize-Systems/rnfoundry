function FRAC=freduce(FRAC)
% freduce: converts a fraction object to lowest common form
%
% FRAC is the fraction (K+N/D) created using fr(...).
% FRAC is reduced if 0<N<D and gcd(N,D)=1
%      or if D=0 and K=0 and N is -1 or 0 or 1 (-inf,nan,inf)
%
% e.g. 0+4/3 -> 1+1/3
%      0-4/3 -> -2+2/3
%
% see also: fr

% Author: Ben Petschel 25/7/2009
%
% Version history:
%   25/7/2009 - first release
%   20/8/2009 - small speedup by using "mygcd"
%   6/10/2012 - replaced 3 instances of "/" by "./"

if ~isa(FRAC,'fr'),
  error('fr:freduce:nonfrac','input argument FRAC is not a fraction object');
end;

for i=1:numel(FRAC),
  K=FRAC(i).whole;
  N=FRAC(i).numer;
  D=FRAC(i).denom;

  % assume K,N,D are all same type
  if isa(K,'double') && max(abs([K,N,D]))>=flintmax
    warning('one of the fraction components is a double exceeding flintmax; results may be inaccurate - use vpi-based fractions instead');
  end;

  s=sign(D);
  fupd=false; % tells if FRAC was updated

  if s==0
    % inf or NaN
    
    % this might cause problems with the class definitions
    K=0;
    N=sign(N);
    D=0;
    fupd=true;
    
  else
    if s==-1
      % bring negative sign to top
      N=-N;
      D=-D;
      fupd=true;
    end;
    if N>=D || N<0
      % reduce N mod D
      R=mod(N,D);
      K=K+(N-R)./D;
      N=R;
      fupd=true;
    end;
    % reduce common factors
    %G=gcd(N,D);
    G=mygcd(N,D);
    if G>1
      N=N./G;
      D=D./G;
      fupd=true;
    end;
  end;

  % update values of FRAC, if necessary
  if fupd
    FRAC(i).numer=N;
    FRAC(i).denom=D;
    FRAC(i).whole=K;
  end;
end; % for i=1:numel(FRAC)

end % main function freduce(...)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function g=mygcd(a,b)
% simple gcd function avoids overhead due to arg checks etc
a=abs(a);
b=abs(b);
while a~=0 && b~=0
  if a>=b
    a=mod(a,b);
  else
    b=mod(b,a);
  end;
end;
if a==0
  g=b;
else
  g=a;
end;

end % helper function mygcd(...)
