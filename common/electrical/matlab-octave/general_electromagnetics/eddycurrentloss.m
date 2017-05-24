function [loss, d] = eddycurrentloss(width, thickness, Bpeak, freq, rho, mu_r, depth)
% eddycurrentloss: calculates the eddy current losses in a flat plate in a
% changing magnetic field. The plate width should be much larger than the
% thickness.
%
%
% Syntax
%
% [loss, d] = eddycurrentloss(width, thickness, Bpeak, freq, rho, mu_r)
% [loss, d] = eddycurrentloss(width, thickness, Bpeak, freq, rho, mu_r, depth)
% 
% Description
% 
% eddycurrentloss calculates the eddy current losses (and skin depth) in a
% flat plate in a changing magnetic field perpendicular to its face. It is
% assumed that the the plate width is much larger than the thickness.
%
%                                                   _
%       /                width                    / .|
%     </---------------------------------------> / . depth
%     /_________________________________________/ .
%     |                    x    x    x          |  ^
%     |                                         |  :
%     |                    x    B    x          |  :
%     |                                         |  : thickness
%     |                    x    x    x          | /:
%     |_________________________________________|/ v
%
%
% The formula used to determine the losses is
%
% depth * width * (thickness./2)^3 * Bpeak^2 * omega^2 / (3*rho)       (1)
%
% taken from [1], where the width and thickness are the plate face
% dimensions and the depth is either the depth of penetration of the
% magnetic field, according to the properties of the material (depth
% calculated using the formula in (2)), or a depth provided by the user. It
% is assumed that the magnetic field is sinusoidal AC field.
% 
% sqrt(2 .* rho / (mu_r .* mu_0 .* 2 .* pi .* freq));                 (2)
%
% [loss, d] = eddycurrentloss(width, thickness, Bpeak, freq, rho, mu_r)
% calculates the skin depth and uses this value in the loss calculation.
% The inputs in htis case are
%
%   width - the width of the plate
%
%   thickness - the thickness of the plate
%
%   Bpeak - the magnitude of the applied sinusoidal AC magnetic field, 
%
%   freq - the frequency in Hertz of the sinusodal variation
%
%   rho - the resistivity of the material
% 
%   mu_r - the relative permeability of the material
%
% [loss, d] = eddycurrentloss(width, thickness, Bpeak, freq, rho, mu_r, depth)
% uses the depth supplied in 'depth'. In this case the value in 'mu_r' is
% ignored and may be empty.
%
%
% [1] G. W. Carter, The Electromagnetic Field in its Engineering Aspects,
% Second Edition. American Elsevier Publishing Company, Inc, 1967.
%

    if nargin == 7
        d = depth;
    elseif nargin == 6
        % first calculate the skin depth
        d = sqrt(2 .* rho / (mu_r .* mu_0 .* 2 .* pi .* freq));
    else
        error('Incorrect number of arguments.')
    end
    
    % now calculate the losses
    loss = d * width .* (thickness./2).^3 .* Bpeak.^2 .* (2 .* pi .* freq).^2 ./ (3 .* rho);

end