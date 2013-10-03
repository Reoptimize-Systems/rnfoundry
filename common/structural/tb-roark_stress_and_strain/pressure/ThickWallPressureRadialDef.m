function delR = ThickWallPressureRadialDef(vars, method)
% function: ThickWallPressureRadialDef
%
% A function for calculating the deformation in a thick-walled pressure
% vessel.
%
% Input:
%
%   vars - (n x p) column vector of values necessary for caluculating the
%          deformation according to the method described in 'method'
%
%   method - string describing the method by which the pressure is to be
%            calculated. These should correspond to the appropriate table 
%            in Roark's Formulas for Stress & Strain.
% 
% Output:
%
%   delR - (n x 2) column vector of values of radial displacements
%          of the circumference of the cylinder, calulated according to
%          'method' the first value is the change in outer radius, the
%          second the change in inner radius
%
    switch method
        
        case '32.1b'
            % Uniform internal pressure q in all directions, ends capped,
            % for a disk or shell
            delR = Table32r1bDelR(vars);
        
        case '32.1c'
            % Uniform external radial pressure q with zero or externally
            % balanced longitudinal pressure
            delR = Table32r1cDelR(vars);
        
        case '32.1f'
            % Linearly varying radial body force from db at centre to zero
            % at outer radius
            delR = Table32r1fDelR(vars);
            
        case '32.1e'
            % Uniform radial body force from db 
            delR = Table32r1eDelR(vars);
            
            
        otherwise
            
            feval(method, vars);
            
    end
    
end