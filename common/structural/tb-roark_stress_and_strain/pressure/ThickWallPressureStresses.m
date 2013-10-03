function sigma = ThickWallPressureStresses(vars, method)
% function: ThickWallPressureStresses
%
% A function for calculating the deformation in a thick-walled pressure
% vessel.
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
        
        case '32.1b'
            % Uniform internal pressure q in all directions, ends capped,
            % for a disk or shell
            sigma = Table32r1bNormStresses(vars);
        
        case '32.1c'
            
            % Uniform external pressure, zero or externally balanced
            % longitudinal pressure
            sigma = Table32r1cNormStresses(vars);
        
        case '32.1e'
            
            % Uniform radial body force from delta
            sigma = Table32r1eNormStresses(vars);
        
        case '32.1f'
            % Linearly varying radial body force from db at centre to zero
            % at outer radius
            sigma = Table32r1fNormStresses(vars);
            
        otherwise
            
            feval(method, vars);
            
    end
    
end