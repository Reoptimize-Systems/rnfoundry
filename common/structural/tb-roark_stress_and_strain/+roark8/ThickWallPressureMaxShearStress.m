function sigma = ThickWallPressureMaxShearStress (vars, method)
% Calculates maximum shear stress in a thick-walled pressure vessel.
%
% Syntax
%
% sigma = roark8.ThickWallPressureMaxShearStress (vars, method)
%
% Input:
%
%   vars - (n x p) column vector of values necessary for caluculating the
%     stresses according to the method described in 'method'
%
%   method - string describing the method by which the pressure is to be
%     calculated. These should correspond to the appropriate table in
%     Roark's Formulas for Stress & Strain.
% 
% Output:
%
%   sigma - (n x 3) column vector of values of sigma, the normal stresses
%     in the longitudinal, circumferential and radial directions
%     respectively.
%

    switch method
        
        case {'13.5.1.1a', 'UniformInternalRadialPress'}
            % Uniform internal radial pressure q, longitudinal pressure
            % zero or externally balanced; for a disk or a shell
            sigma = roark.PressureVessels.ThickWalled.Cylindrical.UniformInternalRadialPress.MaxShearStress (vars);
            
        case {'13.5.1.1c', 'UniformExternalPressLongZero'}
            % Uniform external pressure, zero or externally balanced
            % longitudinal pressure
            sigma = Table32r1cMaxShearStress(vars);
            
        case {'13.5.1.1e'}
            % Uniform radial body force from delta
            sigma = Table32r1eMaxShearStress(vars);
        
        case {'13.5.1.1f'}
            % Linearly varying radial body force from d_b at centre to zero
            % at outer radius
            sigma = Table32r1fMaxShearStress(vars);
            
        otherwise
            
            feval(method, vars);
            
    end
    
end