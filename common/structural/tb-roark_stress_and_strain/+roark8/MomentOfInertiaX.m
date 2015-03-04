function I = MomentOfInertiaX(IVars, IMethod)
% calculates the second moment of inertia of a beam about the x axis.
%
% Syntax
%
% I = MomentOfInertiaX (IVars, IMethod)
%
% Input:
%
%   IVars - (n x p) column vector of values necessary for calculating the
%     second moment of inertia according to the method described in
%     'method'. See the appropriate function for details of the required
%     variables, and exact format of vars.
%
%   IMethod - string describing the method by which the second moment of 
%     inertia is to be calculated. These should correspond to the
%     appropriate table in Roark's Formulas for Stress & Strain, 8th
%     edition. If not one of the listed functions, MomentOfArea will
%     interpret this as a function name and attmpt to evaluate the named
%     function using IVars as the input parameters. Aliases with human
%     understandable names are also available for most sections. Currently,
%     the following sections are available:
%
%           method string options             IVars
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
% Output:
%
%   I - (n x 1) column vector of values of I, the second moment of area of 
%     the cross-section about the x axis, calulated according to 'IMethod'
%
%
% See also: roark8.MomentOfInertiaY
%

    switch IMethod
        
        case {'A1.2', 'Rectangular'}
            % Rectangular
            I = Table1r2I_1(IVars);
            
        case {'A1.3', 'BoxSection'}
            % Hollow Rectangle
            I = Table1r3I_1(IVars);
        
        case {'A1.6', 'IBeam'}
            % I-Beam
            I = Table1r6I_1(IVars);

        case {'A1.15', 'SolidCircle'}
            % Solid Circle
            I = Table1r14I_1(IVars);

        case {'A1.16', 'Annulus'}
            % Hollow Circle
            I = Table1r15I_1(IVars);

        case {'A1.17', 'ThinAnnulus'}
            % Hollow Circle with very thin annulus
            I = Table1r16I_1(IVars);
            
        case {'A1.21', 'HollowCircleSector', 'AnnularSector'}
            % sector of a hollow circle
            I = Table1r21I_1(IVars);

        case 'polygon'

            % use polygeom to calculate the area properties
            for ind = 1:size(IVars, 3)
                [ ~, INER, ~ ] = polygeom( IVars(:,1,ind), IVars(:,2,ind) );
                I(ind,1) = INER(4);
            end

        otherwise

            I = feval(IMethod, IVars);
            
    end

end