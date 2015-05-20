function x_c = DistanceToExtremetiesX (IVars, IMethod)
% calculates the maximum distance from the neutral axis to the extremeties
% of a section in the X direction (i.e. from neutral axis parallel to y
% axis).
%
% Syntax
%
% x_c = DistanceToExtremetiesX (IVars, IMethod)
%
% Input:
%
%   IVars - (n x p) column vector of values necessary for calculating the
%     second moment of inertia according to the method described in
%     'method'. See the appropriate function for details of the required
%     variables, and exact format of vars.
%
%   IMethod - string describing the method by which the distance to the
%     extremities is to be calculated. These should correspond to the
%     appropriate table in Roark's Formulas for Stress & Strain, 8th
%     edition. If not one of the listed functions, MomentOfArea will
%     interpret this as a function name and attmpt to evaluate the named
%     function using IVars as the input parameters. Aliases with human
%     understandable names are also available for most sections. Currently,
%     the following sections are available:
%
%           method string options             IVars
%
%     'A1.15', 'SolidCircle','Circle'      [ R ]
%     'A1.16', 'Annulus', 'HollowCircle'   [ R, Ri ]
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
%   x_c - (n x 1) column vector of values of x_c, the distance to the
%     extremeties in the x direction, calculated according to 'IMethod'
%
%
% See also: roark8.MomentOfInertiaY
%

    switch IMethod
        
%         case {'A1.2', 'Rectangular'}
%             % Rectangular
%             y_c = Table1r2I_1(IVars);
%             
%         case {'A1.3', 'BoxSection', 'HollowRectangle'}
%             % Hollow Rectangle
%             y_c = Table1r3I_1(IVars);
%         
%         case {'A1.6', 'IBeam'}
%             % I-Beam
%             y_c = Table1r6I_1(IVars);

        case {'A1.15', 'SolidCircle', 'Circle'}
            % Solid Circle
            x_c = roark.Sections.DistToExtremities.CircleX (IVars);

        case {'A1.16', 'Annulus', 'HollowCircle'}
            % Hollow Circle
            x_c = roark.Sections.DistToExtremities.HollowCircleX (IVars);

%         case {'A1.17', 'ThinAnnulus'}
%             % Hollow Circle with very thin annulus
%             y_c = Table1r16I_1(IVars);
%             
%         case {'A1.21', 'HollowCircleSector', 'AnnularSector'}
%             % sector of a hollow circle
%             y_c = Table1r21I_1(IVars);
% 
        case 'polygon'

            % use polygeom to calculate the area properties
            for ind = 1:size(IVars, 3)
                thisvars = IVars ( (~isnan (IVars(:,1,ind)) & ~isnan (IVars(:,1,ind))),:,ind);
                [ GEOM, ~, ~ ] = polygeom( thisvars(:,1), thisvars(:,2) );
                % get the distance from the neutral axis by finding the max
                % distance from all y coordinates
                x_c(ind,1) = max ( abs (thisvars(:,1) - GEOM(2)) );
            end

        otherwise

            x_c = feval(IMethod, IVars);
            
    end

end