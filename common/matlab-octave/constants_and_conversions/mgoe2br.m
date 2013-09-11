function Br = mgoe2br(MGOe)
% converts a permanent magnets energy product spec into the equivalent
% remanent flux density
%
% Input
%
%   MGOe - matrix of energy product values to be converted into equivalent
%          magnet coercivity values
%
% Output
%
%   Br - matrix of remanance values for given MGOe values

% Created by Richard Crozier 

    Br = 5e5 .* sqrt(MGOe) .* mu_0 ./ pi ;
    
end