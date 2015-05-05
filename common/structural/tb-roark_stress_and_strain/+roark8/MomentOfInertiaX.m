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
%     'A1.2', 'Rectangular'             [ d, b ]
%     'A1.3', 'BoxSection',             
%           'HollowRectangle'           [ d, b, di, bi ]
%     'A1.6', 'IBeam'                   [ b, t, tw, d ]
%     'A1.15', 'SolidCircle','Circle'   [ R ]
%     'A1.16', 'Annulus'                [ R, Ri ]
%     'A1.17', 'ThinAnnulus'            [ R, t ]
%     'A1.21', 'HollowCircleSector',
%           'AnnularSector'             [ R, t, alpha ]
%     'polygon'                         n x 2 x p matrix 
%
%     The 'polygon' options allows you to supply one or more sets of
%     coordinates representing the border of a closed polygon. These must
%     be the points of the polygon going clockwise around the shape. In the
%     simplest case, a single (n x 2) matrix may be supplied, where each
%     row is an x and y coordinate. Multiple sets of polygons can be
%     supplied using the third dimension of the matrix. Polygons with
%     different numbers of coordinates can be used by filling with NaN
%     values. Rows in each coordinate matrix containing Nan values are
%     ignored.
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
            
        case {'A1.3', 'BoxSection', 'HollowRectangle'}
            % Hollow Rectangle
            I = Table1r3I_1(IVars);
        
        case {'A1.6', 'IBeam'}
            % I-Beam
            I = Table1r6I_1(IVars);

        case {'A1.15', 'SolidCircle', 'Circle'}
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
                thisvars = IVars ( (~isnan (IVars(:,1,ind)) & ~isnan (IVars(:,1,ind))),:,ind);
                [ ~, INER, ~ ] = polygeom( thisvars(:,1), thisvars(:,2) );
                I(ind,1) = INER(4);
            end

        otherwise

            I = feval(IMethod, IVars);
            
    end

end