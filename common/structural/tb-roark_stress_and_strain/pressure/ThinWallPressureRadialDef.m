
function delR = ThinWallPressureRadialDef(vars, method)
% function: ThinWallPressureRadialDef
%
% A function for calculating the deformation in a thin-walled pressure
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
%   delR - (n x 1) column vector of values of delR, the radial displacement
%   of the circumference of the cylinder, calulated according to 'method'
%
    switch method
        
        case '28r1a'
            
            delR = Table28r1aDelR(vars);
            
        case '28r1b'
            
            delR = Table28r1bDelR(vars);
            
        otherwise
            
            error(['method: ', method,', not yet included'])
            
    end
    
end