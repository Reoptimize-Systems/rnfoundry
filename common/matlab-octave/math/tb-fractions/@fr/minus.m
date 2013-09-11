function FRAC = minus(FRAC1,FRAC2)
% fr/minus: Subtracts two fraction objects, or one fraction and a constant
% usage: FRAC = FRAC1 - FRAC2;
% 
% arguments:
%  FRAC,FRAC1,FRAC2 - fraction objects or scalars
%                can be arrays provided the sizes are compatible
%
%
%  See also: uminus, plus

% Author: Ben Petschel 25/7/09
%
% Version history:
%   25/7/09 - first release (using vpi/minus as a template)



% probably should implement fully rather than using plus and uminus, since
% it is used by the comparison operations
FRAC = FRAC1 + uminus(FRAC2);
