function J = PolarMomentOfInertia (vars, JMethod)
% calculates the polar moment of inertia of a section
%
% Syntax
%
% J = roark.PolarMomentOfInertia (vars, JMethod)
%
% Input:
%
%   vars - (n x p) column vector of values necessary for calculating the
%     second moment of inertia according to the method described in
%     'method'. See the appropriate function for details of the required
%     variables, and exact format of vars.
%
%   JMethod - string describing the method by which the second moment of 
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
%     'A1.15', 'SolidCircle','Circle'   [ R ]
%     'A1.16', 'Annulus'                [ R, Ri ]
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
%   J - (n x 1) column vector of values of I, the second moment of area of 
%     the cross-section about the x axis, calulated according to 'IMethod'
%
%
% See also: roark.PolarMomentOfInertia
%

    switch JMethod

        case {'10.1.1', 'SolidCircle', 'Circle'}
            % Solid Circle
            J = roark.PolarMoment.SolidCircle (vars);

        case {'10.1.10', 'Annulus', 'HollowCircle'}
            % Hollow Circle
            J = roark.PolarMoment.HollowCircle (vars);

        case 'polygon'

            % use polygeom to calculate the area properties
            for ind = 1:size(vars, 3)
                thisvars = vars ( (~isnan (vars(:,1,ind)) & ~isnan (vars(:,1,ind))),:,ind);
                [ ~, ~, CPMO ] = polygeom( thisvars(:,1), thisvars(:,2) );
                J(ind,1) = CPMO(5);
            end

        otherwise

            J = feval(IMethod, IVars);
            
    end

end