function [Def, I] = BeamDeflectionY(Ivars, Yvars, E, x, IMethod, beamMethod)
% calulates the deflection of a beam. in the Y axis
%
% Syntax
%
% Def = BeamDeflectionY (Ivars, Yvars, E, x, IMethod, beamMethod)
%
% Input:
%
%   Ivars - (n x p) matrix of values necessary for caluculating the
%     second moment of inertia according to the method described in
%     'IMethod'. See the appropriate function for details of the required
%     variables, and exact format of Ivars.
%
%   Yvars - (n x p) matrix of values necessary for caluculating the
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
%     appropriate table in Roark's Formulas for Stress & Strain, 8th
%     edition. If not one of the listed functions, MomentOfInertiaY will
%     interpret this as a function name and attmpt to evaluate the named
%     function using IVars as the input parameters. Aliases with human
%     understandable names are also available for most sections. Currently,
%     the following sections are available:
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
%     in Roark's Formulas for Stress & Strain, 8h edition.
%
% Output:
%
%   Def - (n x 1) column vector of values of the deflection of a beam, at
%     the positions specified in 'x', calulated according to 'IMethod' and
%     'beamMethod'.
%
%   I - moment of inertia calculated for the beam
%
%

    I = roark8.MomentOfInertiaY (Ivars, IMethod);

    switch beamMethod

        case {'8.1.1f', 'LGRSP', 'lgrsp'}
            % R6 T3.1f
            % Left end guided, right end simply supported (Cantilever),
            % point load
            Def = Table3r1fDef(Yvars, E, I, x);

        case {'8.1.2d', 'LFRFD', 'lfrfd'}
            % R6 T3.2a
            % Left end fixed, right end fixed, distributed force
            Def = Table3r2dDef(Yvars, E, I, x);
            
        case {'8.1.2e', 'LSRSD', 'lsrsd'}
            % R6 T3.2e
            % Left end Simply Supported, right end simply supported,
            % distributed force
            Def = Table3r2eDef(Yvars, E, I, x);
            
        case {'8.1.4d', 'LFRFA', 'lfrfa'}
            % R6 T3.4d
            % Left end fixed, right end fixed, externally created angular
            % displacement
            Def = Table3r4dDef(Yvars, E, I, x);
            
        case {'8.8.2e', 'LSRSDAL', 'lsrsdal'}
            % R6 T10.2e
            % Left end Simply Supported, right end simply supported,
            % distributed force and axial load
            Def = Table10r2eDef(Yvars, E, I, x);
            
        case {'8.8.2f', 'LGRSDAL', 'lgrsdal'}
            % R6 T10.2f
            % Left end guided, right end simply supported,
            % distributed force and axial load
            Def = Table10r2fDef(Yvars, E, I, x);
            

        otherwise
            
            Def = feval(beamMethod, Yvars, E, I, x);

    end
    
end