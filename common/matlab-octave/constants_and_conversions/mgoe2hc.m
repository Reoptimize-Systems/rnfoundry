function Hc = mgoe2hc(MGOe)
% converts a permanent magnets energy product spec into the equivalent
% coercivity
%
% Input
%
%   MGOe - matrix of energy product values to be converted into equivalent
%          magnet coercivity values
%
% Output
%
%   Hc - matrix of coercivity values for given MGOe values

% Created by Richard Crozier 

    Hc = 5e5 .* sqrt(MGOe) ./ pi;
    
end