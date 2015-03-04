function I = MomentOfInertiaY(IVars, IMethod)
% calcualtes the second moment of inertia of a section about the Y axis.
%
% Syntax
%
% I = MomentOfInertiaY(IVars, IMethod)
%
% Input:
%
%   IVars - (n x p) column vector of values necessary for calculating the
%     second moment of inertia according to the method described in 
%     'method'. See the appropriate function for details of the
%     required variables, and exact format of vars.
%
%   IMethod - string describing the method by which the second moment of 
%     inertia is to be calculated. These should correspond to the
%     appropriate table in Roark's Formulas for Stress & Strain, 8th
%     edition. If not one of the listed functions, MomentOfArea will
%     interpret this as a function name and attempt to evaluate the named
%     function using Ivars as the input parameters. Aliases with human
%     understandable names are also available for most sections. Currently,
%     the following sections are available:
%
%           method string                     IVars
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
% Output:
%
%   I - (n x 1) column vector of values of I, the second moment of area of 
%     the cross-section about the y axis, calulated according to 'IMethod'
%
%
% See also: roark8.MomentOfInertiaX
%

    switch IMethod
        
        case {'A1.2', 'Rectangular'}
            % Rectangular
            I = Table1r2I_2(IVars);
            
        case {'A1.3', 'BoxSection'}
            % Hollow Rectangular
            I = Table1r3I_2(IVars);

        case {'A1.6', 'IBeam'}
            % I-Beam
            I = Table1r6I_2(IVars);

        case {'A1.15', 'SolidCircle'}
            % Solid Circle: same as I_1
            I = Table1r14I_2(IVars);

        case {'A1.16', 'Annulus'}
            % Hollow Circle: same as I_1
            I = Table1r15I_2(IVars);

        case {'A1.17', 'ThinAnnulus'}
            % Hollow Circle with very thin annulus: same
            % as I-1
            I = Table1r16I_2(IVars);

        otherwise

            I = feval(IMethod, IVars);
            
    end

end