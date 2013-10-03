function Def = CircPlateDeflection(Gvars, Lvars, E, v, r, plateMethod)
% function for calculating the deflection of a plate
%
% Syntax
%
% Def = CircPlateDeflection(Gvars, Lvars, E, v, r, plateMethod)
%
% Input:
%
%   Gvars - (n x p) matrix of values describing the geometry of the
%            plate(s)
%
%   Lvars - (n x p) matrix of load values necessary for calculating the
%          plate deflection according to the method described in 
%          'plateMethod'. See the appropriate function for details of the
%          required variables, and exact format of Lvars.
%
%   E - Young's modulus of the plate material
%
%   v - poissons ratio for the plate material
%
%   r - row vector of radial position values at which the deflection is to
%       be calculated
%
%   plateMethod - string describing the method by which the plate
%                deflection is to be calculated. These should correspond to
%                the appropriate table in Roark's Formulas for Stress &
%                Strain.
%
% Output:
%
%   Def - (n x 1) column vector of values of the deflection of a dircular plate, at
%         the positions specified in 'x', calulated according to 'IMethod'
%          and 'beamMethod'.
%            
    
    switch plateMethod
        
        case '11.2l'
            % Outer edge free, inner edge fixed
            Def = Table11r2lDef(Gvars, Lvars, E, v, r);
            
        otherwise
            
            Def = feval(plateMethod, Gvars, Lvars, E, v, r);

    end
    
end