function [histloss, eddyloss, excessloss] = ...
    softferrolossrectregionpartcalc(Bx, By, Bz, Hx, Hy, Hz, kc, ke, beta, dx, dy, dz)
% part calculates the values of the hysteresis, eddy current and excess
% loss in a cuboidal region of soft ferromagnetic material
%
% Syntax
%
% [histloss, eddyloss, excessloss] = ...
%       softferrolossrectregpartcalc(Bx, By, Bz, Hx, Hy, Hz, kc, ke, beta, dx, dV)
%
% [histloss, eddyloss, excessloss] = ...
%       softferrolossrectregpartcalc(Bx, By, Bz, Hx, Hy, Hz, kc, ke, beta, dx, dy, dz)
%
% Input
%
%  Bx - matrix of average values of the x directed component of B in a
%   volume of material. The volume is defined by the contents of dx, dy and
%   dz (see their descriptions for details). 
% 
%  By - matrix of average values of the y directed component of B in a
%   volume of material. The volume is defined by the contents of dx, dy and
%   dz (see their descriptions for details). 
% 
%  Bz -  matrix of average values of the z directed component of B in a
%   volume of material. The volume is defined by the contents of dx, dy and
%   dz (see their descriptions for details). 
% 
%  Hx - matrix of average values of the x directed component of H in a
%   volume of material. The volume is defined by the contents of dx, dy and
%   dz (see their descriptions for details).  
% 
%  Hy - matrix of average values of the y directed component of H in a
%   volume of material. The volume is defined by the contents of dx, dy and
%   dz (see their descriptions for details). 
% 
%  Hz - matrix of average values of the z directed component of H in a
%   volume of material. The volume is defined by the contents of dx, dy and
%   dz (see their descriptions for details). 
%
%  kc - 
% 
%  ke -  
% 
%  beta - 
%
%  dx - a scalar value of the spacing between elements in the x direction
% 
%  dy -  
% 
%  dz - 
%
% Output
%
%  histloss -  
% 
%  eddyloss -  
% 
%  excessloss - 
%
% Description
%
% Calculates part of the formula required to perform loss calculations
% using the time-domain method laid out in [1]. The resulting formula
% assumes the direction of motion is constrained so that dB/dt can be
% reasonably represented as:
%
% dB/dx * dx/dt = dB/dx * velocity
%
% [1] D. Lin, P. Zhou, W. N. Fu, Z. Badics, and Z. J. Cendes, A Dynamic
% Core Loss Model for Soft Ferromagnetic and Power Ferrite Materials in
% Transient Finite Element Analysis, IEEE Transactions on Magnetics, vol.
% 40, no. 2, pp. 1318--1321, Mar. 2004.
%
% 

% Copyright Richard Crozier 2012-2013

    if nargin == 11
        dV = dy; % 11th arg is actually volume of elements
    elseif nargin == 12
        dV = dx .* dy .* dz; % dimensions of elements provided
    else
        error('Incorrect number of arguments.')
    end
    
    % the direction of motion is assumed to be in the x direction.
    % Therefore we first the derivatives of each B component with respect
    % to this direction. dx is assumed to be a scalar, therefore Bx, By and
    % Bz must all be sampled at
    dBxVdx = gradient(Bx, dx);
    dByVdx = gradient(By, dx);
    dBzVdx = gradient(Bz, dx);
    
    % generate the part calculation of the hysteresis losses for the
    % region, these must be multiplied by the velocity to get the actual
    % losses
    histloss = abs(Hx .* dBxVdx).^(2/beta) ...
               + abs(Hy .* dByVdx).^(2/beta) ...
               + abs(Hz .* dBzVdx).^(2/beta);
        
    histloss = dV .* (histloss .^ (beta / 2));
    
    histloss = sum(histloss(:));
    
    % generate the part calculation of the eddy current losses for the
    % region, these must be multiplied by the square of the velocity (v^2)
    % to get the actual losses
    Bderivsquares = realpow(dBxVdx, 2) + realpow(dByVdx, 2) + realpow(dBzVdx, 2);
    
    eddyloss =  dV .* (kc / (2 * pi^2)) .* Bderivsquares;
    
    eddyloss = sum(eddyloss(:));
    
    % generate the part calculation of the excess (sometimes innacurately
    % known as anomalous) losses for the region, these must be multiplied
    % by the cube of the square root of the velocity (sqrt(v)^3) to get the
    % actual losses
    
    % Ce is a constant found by the numerical solution of 
    % <latex>
    % \[ Ce = (2*\pi)^{1.5} \frac{2}{\pi} \int_{0}^{\frac{\pi}{2}} \cos^{1.5}(\theta) d \theta \]
    % <latex>
    % This is described in [1]
    Ce = 8.763363;
    
    excessloss = dV .* (ke / Ce) .*  realpow(Bderivsquares, 0.75);
    
    excessloss = sum(excessloss(:));

end


