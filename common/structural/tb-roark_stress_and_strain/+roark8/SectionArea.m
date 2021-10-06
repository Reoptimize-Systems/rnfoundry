function A = SectionArea (vars, method)
% calculates the cross-sectional area of a section as tabulated in Roark's
% Formulas for Stress and Strain 8th Edition
%
% Syntax
%
% A = SectionArea (vars, method)
%
% Input
%
%   vars - (n x p) column vector of values necessary for calculating the
%     section area according to the method described in 'method'. See the
%     appropriate function for details of the required variables, and exact
%     format of vars.
%
%   method - string describing the method by which the cross-sectional area
%     is to be calculated. These should correspond to the appropriate table
%     in Roark's Formulas for Stress & Strain, 8th edition. If not one of
%     the listed functions, SectionArea will interpret this as a function
%     name and attmpt to evaluate the named function using IVars as the
%     input parameters. Aliases with human understandable names are also
%     available for most sections. Currently, the following sections are
%     available:
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
%   A - (n x 1) column vector of values of A, the cross-sectional area of
%     the section, calulated according to 'method'
%

    switch method
        
        case {'A1.2', 'Rectangular'}
            % Rectangular
            A = Table1r2Area(vars);
            
        case {'A1.3', 'BoxSection', 'HollowRectangle'}
            % Hollow Rectangle
            A = Table1r3Area(vars);
        
        case {'A1.6', 'IBeam'}
            % I-Beam
            A = Table1r6Area(vars);

        case {'A1.8', 'UnequalAngle', 'LProfile'}
            % Unequal-Legged Angle (L profile)
            A = roark.Sections.Areas.UnequalLegAngle (vars);
            
        case {'A1.15', 'SolidCircle', 'Circle'}
            % Solid Circle
            A = Table1r14Area(vars);

        case {'A1.16', 'Annulus'}
            % Hollow Circle
            A = Table1r15Area(vars);

        case {'A1.17', 'ThinAnnulus'}
            % Hollow Circle with very thin annulus
            A = Table1r16Area(vars);
            
%         case {'A1.21', 'HollowCircleSector', 'AnnularSector'}
%             % sector of a hollow circle
%             I = Table1r21I_1(IVars);

        case 'polygon'

            % use polygeom to calculate the area properties
            for ind = 1:size(vars, 3)
                thisvars = vars ( (~isnan (vars(:,1,ind)) & ~isnan (vars(:,1,ind))),:,ind);
                [ GEOM, ~, ~ ] = polygeom(thisvars(:,1), thisvars(:,2));
                A(ind,1) = GEOM(1);
            end
            
        otherwise

            A = feval(method, vars);
            
    end
end