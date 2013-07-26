function d = skindepth(rho, mu_r, freq)
% calculates the skin depth at a given frequency of oscillating magnetic
% field in a good conductor
%
% Syntax
%
% d = skindepth(rho, mu_r, freq)
%
% Input
%
% rho - the resistivity of the material 
%
% mu_r - relative permeability of the material
%
% freq - frequency of field oscillation in Hz
%
% Output
%
% d - the effective depth of penetration of the the field (the depth where
%   its value has fallen to around 1/e of that on the surface)
%

    d = sqrt(2 .* rho ./ (mu_r .* mu_0 .* 2 .* pi .* freq));
    
    
end