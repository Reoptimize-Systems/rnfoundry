function S = sign(FRAC)
% fr/sign: sign of a fraction object or array
% usage: S = sign(FRAC);
% 
% arguments:
%  FRAC - a fraction object or array
%
%  S   - Takes the sign of FRAC. One of
%        {-1, 0, +1, NaN} will be returned.
%        S will be a double.
%
%  Example:
%  f=fr(0,[-1,0,1],0)  % fractional rep of [-inf,nan,inf]
%  sign(f)             % returns [-1,NaN,1]
%
%  See also: abs

% Author: Ben Petschel 25/7/09
%
% Version history:
%   25/7/09 - first release (using vpi/sign as a template)

S = zeros(size(FRAC));

fnnan = ~isnan(FRAC); % boolean true for fractions not NaN
S(~fnnan) = nan; % sign is nan for fractions not not NaN (i.e. NaN)

% sign of finite and inf is sign of the whole part
% or sign of the numerator if the whole part is zero
SK = sign([FRAC(fnnan).whole]);
S(fnnan) = SK;
S(fnnan & (S==0)) = sign([FRAC(fnnan & (S==0)).numer]);

end
