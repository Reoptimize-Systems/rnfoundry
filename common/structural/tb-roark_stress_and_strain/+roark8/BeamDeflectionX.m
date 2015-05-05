function [Def, I, thetaA, MA, RA, yA] = BeamDeflectionX(Ivars, Yvars, E, x, IMethod, beamMethod)
% calculates the deflection of a beam about the x-axis using formulas from
% 'Roark's Formulas for stress and strain'
%
% Syntax 
%
% Def = roark8.BeamDeflectionX (Ivars, Yvars, E, x, IMethod, beamMethod)
%
% Input
%
%   Ivars - (n x p) matrix of values necessary for calculating the
%     second moment of inertia according to the method described in
%     'IMethod'. See the appropriate function for details of the required
%     variables, and exact format of Ivars.
%
%   Yvars - (n x p) matrix of values necessary for calculating the
%     beam deflection according to the method described in 'beamMethod'.
%     See the appropriate function for details of the required variables,
%     and exact format of Yvars.
%
%   E - Young's modulus of the beam material
%
%   x - row vector of position values at which the deflection is to be
%     calculated 
%
%   IMethod - string describing the method by which the second moment of 
%     inertia is to be calculated. These should correspond to the
%     appropriate table in Roark's Formulas for Stress & Strain. If not one
%     of the listed functions, MomentOfInertiaX will interpret this as a
%     function name and attmpt to evaluate the named function using IVars
%     as the input parameters. Aliases with human understandable names are
%     also available for most sections. Currently, the following sections
%     are available:
%
%          method string options            IVars
%
%     'A1.2', 'Rectangular'           [ b, d ]
%     'A1.3', 'BoxSection'            [ b, d, di, bi ]
%     'A1.6', 'IBeam'                 [ b, t, tw, d ]
%     'A1.15', 'SolidCircle'          [ R ]
%     'A1.16', 'Annulus'              [ R, Ri ]
%     'A1.17', 'ThinAnnulus'          [ R, t ]
%     'A1.21', 'HollowCircleSector',
%           'AnnularSector'           [ R, t, alpha ]
%     'polygon'                        n x 2 x p matrix 
%
%   beamMethod - string describing the method by which the beam deflection
%     is to be calculated. These should correspond to the appropriate table
%     in Roark's Formulas for Stress & Strain.
%
% Output
%
%   Def - (n x 1) column vector of values of the deflection of a beam, at
%     the positions specified in 'x', calulated according to 'IMethod' and
%     'beamMethod'.
%            
%   I - moment of inertia calculated for the beam
%
%   thetaA, MA, RA, yA
%
% See also: roark8.MomentOfInertiaX
%

    I = roark8.MomentOfInertiaX (Ivars, IMethod);
    
    switch beamMethod
        
        case {'8.1.1d', 'LFRFP', 'lfrfp'}
            % R6 T3.1d
            % Left end fixed, right end fixed, point load
            Def = Table3r1dDef (Yvars, E, I, x);
        
        case {'8.1.1e', 'LSRSP', 'lsrsp'}
            % R6 T3.1e
            % Left end Simply Supported, right end simply supported,
            % point load
            Def = Table3r1eDef (Yvars, E, I, x);
        
        case {'8.1.1f', 'LGRSP', 'lgrsp'}
            % R6 T3.1f
            % Left end guided, right end simply supported,
            % point load
            Def = Table3r1fDef (Yvars, E, I, x);
        
        case {'8.1.2a', 'LURFD', 'lurfd'}
            % R6 T3.2a
            % Left end free, right end fixed, distributed force
            Def = Table3r2aDef (Yvars, E, I, x);
            
        case {'8.1.2d', 'LFRFD', 'lfrfd'}
            % R6 T3.2a
            % Left end fixed, right end fixed, distributed force
            Def = Table3r2dDef (Yvars, E, I, x);
            
        case {'8.1.2e', 'LSRSD', 'lsrsd'}
            % R6 T3.2e
            % Left end Simply Supported, right end simply supported,
            % distributed force
            [Def, thetaA, MA, RA, yA] = Table3r2eDef (Yvars, E, I, x);
            
        case {'8.1.4d', 'LFRFA', 'lfrfa'}
            % R6 T3.4d
            % Left end fixed, right end fixed, externally created angular
            % displacement
            Def = Table3r4dDef (Yvars, E, I, x);
            
        case {'8.8.2e', 'LSRSDAL', 'lsrsdal'}
            % R6 T10.2e
            % Left end Simply Supported, right end simply supported,
            % distributed force and axial load
            Def = Table10r2eDef (Yvars, E, I, x);
            
        case {'8.8.2f', 'LGRSDAL', 'lgrsdal'}
            % R6 10.2f
            % Left end guided, right end simply supported,
            % distributed force and axial load
            Def = Table10r2fDef (Yvars, E, I, x);
            
        otherwise
            
            Def = feval (beamMethod, Yvars, E, I, x);

    end
    
end