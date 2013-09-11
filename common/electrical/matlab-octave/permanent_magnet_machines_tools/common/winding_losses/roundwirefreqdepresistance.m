function Rac = roundwirefreqdepresistance(a, Rdc, rho, mu_r, freq)
% calcuates the AC resitance of a wire of round cross-section due to the
% skin effect
%
% Syntax
%
% Rac = roundwirefreqdepresistance(a, Rdc, rho, mu_r, freq)
%
% Input
%
%   a - the conductor radius
%
%   Rdc - the DC resistance of the wire
% 
%   rho - the resistivity of the wire material
% 
%   mu_r - the relative permeability of the wire material
% 
%   freq - the frequency of the current waveform
%
% Output
%
%   Rac - the AC wire resistance including the skin effect
%
% Description
%
% The AC winding resistance is calculated according to the formulas
% presented in 'The Analysis of Eddy Currents', Richard L Stoll, Chapter 2,
% Section 2.8, page 25
%
% 

% Copyright Richard Crozier 2012 - 2012

    % determine the skin depth
    delta = skindepth(rho, mu_r, freq);
    
    % calculate the AC resistance, this is dependent on the ratio of the
    % wire radius to the skin depth, as described in 'The Analysis of Eddy
    % Currents', Richard L Stoll, Chapter 2, Section 2.8, page 25
    if a > 7 * delta
        
        Rac = Rdc .* ( (a ./ (2.*delta)) + 0.25 + (3.*delta ./ (32.*a)) );
        
    else
        
        Rac = Rdc .* ( 1 + a.^2 ./ (4.*delta.^2) );
        
    end


end