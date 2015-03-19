function SectionPlot(vars, method)
% calculates the second moment of inertia of a beam about the x axis.
%
% Syntax
%
% SectionPlot (IVars, IMethod)
%
% Input:
%
%   vars - (n x p) column vector of values necessary for calculating the
%     second moment of inertia according to the method described in
%     'method'. See the appropriate function for details of the required
%     variables, and exact format of vars.
%
%   method - string describing the method by which the second moment of 
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
%     The 'polygon' options allows you to supply a set of coordinates
%     representing the border of a closed polygon. These must be the points
%     of the polygon going clockwise around the shape. A single (n x 2)
%     matrix may be supplied, where each row is an x and y coordinate.
%
% Output:
%
%   None
%
% See also: roark8.MomentOfInertiaY
%

    figure;
    
    switch method
        
        case {'A1.2', 'Rectangular'}
            % Rectangular
            figure;
            rectangle('Position',[0,0,vars(2),vars(1)]) ;
            axis equal
            
        case {'A1.3', 'BoxSection', 'HollowRectangle'}
            % Hollow Rectangle
            figure;
            rectangle('Position',[0,0,vars(2),vars(1)]) ; 
            hold on; 
            rectangle('Position',[(vars(2)-vars(4))/2,(vars(1)-vars(3))/2,vars(4),vars(3)]); 
            hold off; 
            axis equal
        
%         case {'A1.6', 'IBeam'}
%             % I-Beam
%             I = Table1r6I_1(IVars);
% 
%         case {'A1.15', 'SolidCircle', 'Circle'}
%             % Solid Circle
%             I = Table1r14I_1(IVars);
% 
%         case {'A1.16', 'Annulus'}
%             % Hollow Circle
%             I = Table1r15I_1(IVars);
% 
%         case {'A1.17', 'ThinAnnulus'}
%             % Hollow Circle with very thin annulus
%             I = Table1r16I_1(IVars);
%             
%         case {'A1.21', 'HollowCircleSector', 'AnnularSector'}
%             % sector of a hollow circle
%             I = Table1r21I_1(IVars);
% 
%         case 'polygon'
% 
%             % use polygeom to calculate the area properties
%             for ind = 1:size(IVars, 3)
%                 thisvars = vars ( (~isnan (vars(:,1,ind)) & ~isnan (vars(:,1,ind))),:,ind);
%                 [ ~, INER, ~ ] = polygeom( thisvars(:,1), thisvars(:,2) );
%                 I(ind,1) = INER(4);
%             end

        otherwise

            error ('Plot for section type not yet implemented.')
            
    end

end