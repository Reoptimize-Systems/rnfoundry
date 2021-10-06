function sigma = ThickWallPressureNormalStresses (vars, method)
% Calculating the normal stresses in a thick walled pressure vessel.
%
% Syntax
%
% sigma = roark8.ThickWallPressureNormalStresses (vars, method)
%
% Input:
%
%   vars - (n x p) column vector of values necessary for calculating the
%     stresses according to the method described in 'method'
%
%   method - string describing the method by which the stress is to be
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
            sigma = roark.PressureVessels.ThickWalled.Cylindrical.UniformInternalRadialPress.NormalStresses (vars);
        
        case {'13.5.1.1b', 'UniformInternalRadialPressEndsCapped'}
            % Uniform internal radial pressure q in all directions, ends
            % capped, for a disk or shell
            sigma = roark.PressureVessels.ThickWalled.Cylindrical.UniformInternalRadialPressEndsCapped.NormalStresses (vars);
        
        case {'13.5.1.1c', 'UniformExternalPressLongZero'}
            % Uniform external radial pressure q, longitudinal pressure
            % zero or externally balanced; for a disk or a shell
            sigma = roark.PressureVessels.ThickWalled.Cylindrical.UniformExternalPressLongZero.NormalStresses (vars);
        
        case {'13.5.1.1e'}
            % Uniform radial body force from delta
            sigma = Table32r1eNormStresses(vars);
        
        case {'13.5.1.1f'}
            % Linearly varying radial body force from db at centre to zero
            % at outer radius
            sigma = Table32r1fNormStresses(vars);
            
        otherwise
            
            if isa (method, 'function_handle') || (ischar (method) && exist (method, 'file'))
                feval(method, vars);
            else
                error ('method not recognised, may not yet be implemented');
            end
            
    end
    
end