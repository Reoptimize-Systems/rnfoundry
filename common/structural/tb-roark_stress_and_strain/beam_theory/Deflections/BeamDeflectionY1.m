function Def = BeamDeflectionY1(Ivars, Yvars, E, x, IMethod, beamMethod)
% function: BeamDeflectionY1
%
% function for calculating the deflection of a beam.
%
% Input:
%
%   Ivars - (n x p) matrix of values necessary for caluculating the
%    second moment of inertia according to the method described in
%    'IMethod'. See the appropriate function for details of the required
%    variables, and exact format of Ivars.
%
%   Yvars - (n x p) matrix of values necessary for calculating the
%    beam deflection according to the method described in 'beamMethod'. See
%    the appropriate function for details of the required variables, and
%    exact format of Yvars.
%
%   E - Young's modulus of the beam material
%
%   x - row vector of position values at which the deflection is to be
%    calculated 
%
%   IMethod - string describing the method by which the second moment of 
%    inertia is to be calculated. These should correspond to the
%    appropriate table in Roark's Formulas for Stress & Strain.
%
%   beamMethod - string describing the method by which the beam deflection
%    is to be calculated. These should correspond to the appropriate table
%    in Roark's Formulas for Stress & Strain. Currently the following
%    strings are possible:
%
%    '3.1d' : Left end fixed, right end fixed, point load
%
%    '3.1e' : Left end Simply Supported, right end simply supported, point 
%       load
%  
%    '3.1f' : Left end guided, right end simply supported, point load
% 
%    '3.2a' : Left end fixed, right end free, distributed force
% 
%    '3.2d' : Left end fixed, right end fixed, distributed force
%
%    '3.2e' : Left end Simply Supported, right end simply supported,
%       distributed force
%
%    '3.4d' : Left end fixed, right end fixed, externally created angular
%       displacement
%  
%    '10.2e' : Left end Simply Supported, right end simply supported,
%       distributed force and axial load
%   
%    '10.2f' : Left end guided, right end simply supported, distributed 
%       force and axial load
%
%    To calculate the deflection, internally, a corresponding function is
%    called, e.g. Table3r1dDef is called when the '3.1' string is used. To
%    understand the correct values and format of Yvars, examine the help
%    for the appropriate fuction.
%
% Output:
%
%   Def - (n x 1) column vector of values of the deflection of a beam, at
%         the positions specified in 'x', calulated according to 'IMethod'
%          and 'beamMethod'.
%            

    I = MomentOfInertiaY1(Ivars, IMethod);
    
    switch beamMethod
        
        
        case '3.1d'
            % Left end fixed, right end fixed, point load
            Def = Table3r1dDef(Yvars, E, I, x);
        
        case '3.1e'
            % Left end Simply Supported, right end simply supported,
            % point load
            Def = Table3r1eDef(Yvars, E, I, x);
        
        case '3.1f'
            % Left end guided, right end simply supported,
            % point load
            Def = Table3r1fDef(Yvars, E, I, x);
        
        case '3.2a'
            % Left end fixed, right end free, distributed force
            Def = Table3r2aDef(Yvars, E, I, x);
            
        case '3.2d'
            % Left end fixed, right end fixed, distributed force
            Def = Table3r2dDef(Yvars, E, I, x);
            
        case '3.2e'
            % Left end Simply Supported, right end simply supported,
            % distributed force
            Def = Table3r2eDef(Yvars, E, I, x);
            
        case '3.4d'
            % Left end fixed, right end fixed, externally created angular
            % displacement
            Def = Table3r4dDef(Yvars, E, I, x);
            
        case '10.2e'
            % Left end Simply Supported, right end simply supported,
            % distributed force and axial load
            Def = Table10r2eDef(Yvars, E, I, x);
            
        case '10.2f'
            % Left end guided, right end simply supported,
            % distributed force and axial load
            Def = Table10r2fDef(Yvars, E, I, x);
            

        otherwise
            
            Def = feval(beamMethod, Yvars, E, I, x);

    end
    
end