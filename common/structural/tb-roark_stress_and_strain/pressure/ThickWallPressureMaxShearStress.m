function sigma = ThickWallPressureMaxShearStress(vars, method)
% ThickWallPressureStresses: calculates the deformation in a thick-walled
% pressure vessel.
%
% Input:
%
%   vars - (n x p) column vector of values necessary for caluculating the
%          stresses according to the method described in 'method'
%
%   method - string describing the method by which the pressure is to be
%            calculated. These should correspond to the appropriate table 
%            in Roark's Formulas for Stress & Strain.
% 
% Output:
%
%   sigma - (n x 3) column vector of values of sigma, the normal stresses in
%           the longitudinal, circumferential and radial directions
%           respectively.
%
    switch method
        
        case '32.1c'
            
            % Uniform external pressure, zero or externally balanced
            % longitudinal pressure
            sigma = Table32r1cMaxShearStress(vars);
            
        case '32.1e'
            
            % Uniform radial body force from delta
            sigma = Table32r1eMaxShearStress(vars);
        
        case '32.1f'
            
            % Linearly varying radial body force from d_b at centre to zero
            % at outer radius
            sigma = Table32r1fMaxShearStress(vars);
            
        otherwise
            
            feval(method, vars);
            
    end
    
end