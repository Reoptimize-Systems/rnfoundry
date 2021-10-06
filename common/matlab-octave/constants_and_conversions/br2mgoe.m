function MGOe = br2mgoe(Br)
% converts a permanent magnets energy product spec into the equivalent
% remanent flux density
%
% Input
%
%   Br - matrix of remanance values for given MGOe values
%
% Output
%
%   MGOe - matrix of energy product values to be converted into equivalent
%          magnet coercivity values
    
% Created by Richard Crozier 

    MGOe = realpow(pi .* Br ./ (5e5 .* mu_0), 2);
    
end