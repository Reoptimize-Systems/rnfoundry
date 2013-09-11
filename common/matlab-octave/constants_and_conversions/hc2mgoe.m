function MGOe = hc2mgoe(Hc)
% converts a permanent magnets coercivity into the equivalent
% energy product
%
% Input
%
%   Hc - matrix of coercivity values for given MGOe values
%
% Output
%
%   MGOe - matrix of energy product values to be converted into equivalent
%          magnet coercivity values
%

% Created by Richard Crozier 

    MGOe = realpow(Hc .* pi ./ 5e5, 2);
    
end