function varargout=rat(FRAC)
% fr/rat - finds the numerator and denominator of a fraction object
%
% usage: [N,D] = rat(FRAC)
% usage: [K,N,D] = rat(FRAC)
% 
%
% arguments: (input)
%  FRAC1,FRAC2 - fraction objects or compatible objects
%      can be arrays if the sizes are compatible
%
% arguments: (output)
%  K,N,D - if called with 3 outputs, returns K,N,D where FRAC=K+N/D
%          if called with 2 outputs, returns N,D where FRAC=N/D
%
% Example:
%  f=fr(1,1,2);   % the fraction 1+1/2 = 3/2
%  [n,d]=rat(f)   % returns n=3, d=2
%  [k,n,d]=rat(f) % returns k=1, n=1, d=2
%
%
%  See also: fr, double

% Author: Ben Petschel 25/7/09
%
% Version history:
%   25/7/09 - first release

if nargout<2 || nargout>3
  error('fr:rat:nargout','must have 2 or 3 output arguments');
end;

if nargin~=1,
  error('fr:rat:nargin','must have exactly 1 input argument');
end;

if numel(FRAC)==1,
  K=FRAC.whole;
  N=FRAC.numer;
  D=FRAC.denom;

else
  K=reshape([FRAC(:).whole],size(FRAC));
  N=reshape([FRAC(:).numer],size(FRAC));
  D=reshape([FRAC(:).denom],size(FRAC));

end;

if nargout==3,
  varargout{1}=K;
  varargout{2}=N;
  varargout{3}=D;
else
  varargout{1}=N+K*D;
  varargout{2}=D;
end;
    
end
