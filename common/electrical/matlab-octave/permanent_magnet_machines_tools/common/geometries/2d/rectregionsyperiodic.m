function [nodes, nodeids, links, rectcentres, spacecentres] = ...
    rectregionsyperiodic(x, y1, y2, xoffset, ydisp, tol, nodecount)
% Generates node locations and segment joints for the drawing of a periodic
% region made up of multiple rectangles
%
% Syntax
%
% [nodes, nodeids, links, rectcentres, spacecentres] = ...
%       rectregionsyperiodic(x, y1, y2, xoffset, ydisp, tol, nodecount)
% 
% Description
%
% The base condition of the rectangles is as shown below. It is assumed
% that there is an infinite series of this pattern in both directions. 
%
%              ***************                                            
%              *             *                                              
%              *             * 
%              ***************                                            
%              *             *                                              
%              *             *                                                                                          
%              *             *                                              
%              *      B      *                                                                                     
%              *             *                                              
%              *             *                                              
%              *             *            
%              ***************   ^                                                                                       
%              *             *   |                                           
%              *             *   | y2
%              *             *   |
%              *             *   v
%              ***************   ^                                           
%              *             *   |
%              *             *   |                                                                                        
%              *             *   | y1                                          
%              *      A      *   |                                     
%              *             *   |               y                           
%              *             *   |              / \
%              *             *   v               |                          
%              ***************                   |                         
%              *             *                   |_____\
%              *             *                         / x
%              ***************                                            
%              <------------->
%                     x
%
% This geometry is drawn in a periodic way in the y direction, 'wrapping'
% around at the top and bottom. The 'A' and 'B' rectangles will always be
% present, but the number of 'spaces' between these may vary depending on
% the tolerance and specified position in the y direction.
%
% Inputs
%
%  x - width of the region in the x direction
%
%  y1 - rectangular region A and B height in y direction
%
%  y2 - height of rectangular region between A and B in y direction
%
%  xoffset - displacement of region center from zero in x direction
%
%  ydisp - displacement in the 'y' direction. The outer boundary of
%    the geometry will not actually change, the gemoetry will be modified
%    as though rotated on a cylindrical surface, wrapping around at the top
%    and bottom. The ydisp units are the same as the y1 and y2 dimensions.
%
%  tol - optional tolerance at which distances bewteen elements are
%    considered to be of size zero. If not supplied, 1e-5 is used.
%
%  nodecount - optional starting count for the node ids used to specify the
%    links between nodes. This number will be added to all the link
%    numbers. Defaults to zero if not supplied.  
%
% Output
%
%  nodes - (n x 2) matrix of x and y coordinates of the geometry
%
%  nodeids -  vector of id numbers for the nodes in nodes
%
%  links - (m x 2) matrix of links bewteen nodes specified as the two id
%    numbers of the nodes
%
%  rectcentres - (p x 2) matrix of x and y coordinates of the centers of
%    the 'A' and 'B' rectangles, or their component parts (can be split
%    into two when 'wrapped')
%
%  spacecentres - (p x 2) matrix of x and y coordinates of the centers of
%    the spaces between the 'A' and 'B' rectangles.
%
%


    if nargin < 6 || isempty(tol)
        tol = 1e-5;
    end

    if nargin < 7
        nodecount = 0;
    end
    
    yperiod = 2*y1 + 2*y2;
    
    % First define the rectangle 1 coordinates, the other rectangels will
    % be defined relative to these
    % Bottom left of bottom rect
	xycoords(1,:) = [-x/2, y2/2];
    % Bottom right of bottom rect
	xycoords(2,:) = [x/2, y2/2];
    % Top left of bottom rect
	xycoords(3,:) = [-x/2, xycoords(1,2) + y1];
    % Top right of bottom rect
    xycoords(4,:) = [x/2, xycoords(2,2) + y1];
    % Top rectangle is bottom rectangle shifted by half period
    xycoords(5:8,1) = xycoords(1:4,1);
    xycoords(5:8,2) = xycoords(1:4,2) + yperiod/2;
    
    % Label Coordinate, bottom rect then top rect
    xycoords(9,:) = [xycoords(1,1)+(xycoords(2,1)-xycoords(1,1))/2, xycoords(1,2)+(xycoords(3,2)-xycoords(1,2))/2];
    xycoords(10,:) = [xycoords(5,1)+(xycoords(6,1)-xycoords(5,1))/2, xycoords(5,2)+(xycoords(7,2)-xycoords(5,2))/2];
    
    % Map coordinates onto cylinder surface
    xycoords(:,2) = xycoords(:,2) .* pi ./ (yperiod/2);
    
    ydisp = ydisp .* pi ./ (yperiod/2);
    
    % rotate around cylinder
    xycoords(:,2) = rem(xycoords(:,2) + ydisp, 2*pi);
    
    % negative coordinates imply they are in fact rotated backwards which
    % we will take to be positions relative to the top of the sim
    xycoords(sign(xycoords(:,2)) == -1, 2) = (2*pi) + xycoords(sign(xycoords(:,2)) == -1, 2); 
     
    % Map back to x-y plane
    xycoords(:,2) = (yperiod/2) .* xycoords(:,2) ./ pi;
    
    minsizes = max(y1 * 0.0005, tol);
    
    bottomnodes = [-x/2, 0; x/2, 0];
    topnodes = [-x/2, yperiod; x/2, yperiod];
    
    % now generate the nodes and joints according to the specified
    % tolerances
    if abs(xycoords(1,2)) <= minsizes || abs(((2*(yperiod/2)) - abs(xycoords(1,2)))) <= minsizes
        % The bottom of the rectangle A is against the base of the sim
        % (or at least within 1e-5 m of it (0.01 mm)
        % Draw line 3--4 and join up nodes to base then draw top mag
        % rectangle 1, link nodes 4-6 and 8-top
        
        %                    9-8
        %            8 *************** 9                                                                                          
        %              *             *                                              
        %          6-8 *             * 7-9
        %              *             *          
        %            6 *************** 7                                            
        %              *     7-6     *                                              
        %              *             *                                              
        %              *             *                                              
        %              *             *                                              
        %          6-4 *      B      * 5-7                                          
        %              *             *                                              
        %              *             *                                              
        %              *             *                                              
        %              *     4-5     *            
        %            4 *************** 5              nodes:                                                                              
        %              *             *                                              
        %          2-4 *             * 3-5      [ bottom left     0
        %              *             *            bottom right    1   
        %            2 *************** 3          rectA top left  2                                                
        %              *     3-2     *            rectA top right 3  
        %              *             *            rectB bot left  4                                                 
        %              *             *            rectB bot right 5                                                
        %              *             *            rectB top left  6                                                
        %          2-0 *      A      * 1-3        rectB top right 7                                             
        %              *             *            top left        8                                                
        %              *             *            top right ]     9                                                
        %              *             *                                                            
        %              *     0-1     *                                                            
        %            0 *************** 1              
        %                    
        
        nodes = [bottomnodes; 
                 xycoords(3:4,:); 
                 xycoords(5:8,:); 
                 topnodes];
             
        links = [0,1; 0,2; 2,3; 1,3; 4,5; 5,7; 7,6; 6,4; 2,4; 3,5; 6,8; 7,9; 8,9]; 

        rectcentres = rectcentre( [nodes(1,:); nodes(5,:)], ...
                                  [nodes(4,:); nodes(8,:)] );
        
        % Add a designation to determine which rectangle parts are which
        rectcentres = [ rectcentres, [0; 1] ];
        
        spacecentres = rectcentre( [nodes(3,:); nodes(7,:)], ...
                                   [nodes(6,:); nodes(10,:)] );

    elseif abs(xycoords(3,2)) <= minsizes || abs(((2*(yperiod/2)) - abs(xycoords(3,2)))) <= minsizes
        % The top of the rectangle A is against the base of the sim
        % (or at least within 1e-5 m of it (0.1 mm), or equivalently, the
        % top of the top rectangle A is against the top of the sim
        % Draw line 1--2 and join up nodes to top then draw top mag
        % rectangle 1, link nodes 8-2 and 6-bot
        
        %            8 *************** 9                                            
        %              *     8-9     *                                              
        %              *             *                                              
        %              *             *                                              
        %              *             *                                              
        %          6-8 *      A      * 7-9                                          
        %              *             *                                              
        %              *             *                                              
        %              *             *                                              
        %              *     6-7     *            
        %            6 *************** 7                                                                                          
        %              *             *                                              
        %          6-4 *             * 7-5
        %              *             *  
        %            4 *************** 5                                            
        %              *     4-5     *
        %              *             *                                              
        %              *             *               nodes:                          
        %              *             *                                              
        %          2-4 *      B      * 3-5      [ bottom left     0                                 
        %              *             *            bottom right    1                                  
        %              *             *            rectB bot left  2                                  
        %              *             *            rectB bot right 3                                  
        %              *     2-3     *            rectB top left  4                                  
        %            2 *************** 3          rectB top right 5                                  
        %              *             *            rectA bot left  6
        %          0-2 *             *  1-3       rectA bot right 7
        %              *             *            top left        8
        %            0 *************** 1          top right ]     9                                  
        %                    0-1

        nodes = [bottomnodes; 
                 xycoords(5:8,:); 
                 xycoords(1:2,:); 
                 topnodes];
             
        links = [0,1; 0,2; 1,3; 2,3; 3,5; 5,4; 4,2; 4,6; 5,7; 6,7; 6,8; 7,9; 8,9]; 

        rectcentres = rectcentre( [nodes(3,:); nodes(7,:)], ...
                                  [nodes(6,:); nodes(10,:)] );
                  
        % Add a designation to determine which rectangle parts are which
        rectcentres = [ rectcentres, [1; 0] ];
        
        spacecentres = rectcentre( [nodes(1,:); nodes(5,:)], ...
                                   [nodes(4,:); nodes(8,:)] );        
        
    elseif abs(xycoords(5,2)) <= minsizes || abs(((2*(yperiod/2)) - abs(xycoords(5,2)))) <= minsizes
        % The bottom of the rectangle B is against the base of the sim
        % (or at least within 1e-4 m of it (0.1 mm)
        % Draw line 7--8 and join up nodes to base then draw bottom mag
        % rectangle 1, link nodes 8-2 and 4-top

        %                    9-8
        %            8 *************** 9                                                                                          
        %              *             *                                              
        %          6-8 *             * 7-9
        %              *             *          
        %            6 *************** 7                                            
        %              *     7-6     *                                              
        %              *             *                                              
        %              *             *                                              
        %              *             *                                              
        %          6-4 *      A      * 5-7                                          
        %              *             *                                              
        %              *             *                                              
        %              *             *                                              
        %              *     4-5     *            
        %            4 *************** 5              nodes:                                                                              
        %              *             *                                              
        %          2-4 *             * 3-5      [ bottom left     0
        %              *             *            bottom right    1   
        %            2 *************** 3          rectB top left  2                                                
        %              *     3-2     *            rectB top right 3  
        %              *             *            rectA bot left  4                                                 
        %              *             *            rectA bot right 5                                                
        %              *             *            rectA top left  6                                                
        %          2-0 *      B      * 1-3        rectA top right 7                                             
        %              *             *            top left        8                                                
        %              *             *            top right ]     9                                                
        %              *             *                                                            
        %              *     0-1     *                                                            
        %            0 *************** 1              
        %            

        nodes = [bottomnodes; 
                 xycoords(7:8,:); 
                 xycoords(1:4,:); 
                 topnodes];
             
        links = [0,1; 1,3; 3,2; 2,0; 3,5; 2,4; 4,5; 5,7; 7,6; 6,4; 7,9; 6,8; 8,9];
        
        rectcentres = rectcentre( [nodes(1,:); nodes(5,:)], ...
                                  [nodes(4,:); nodes(8,:)] );
              
        % Add a designation to determine which rectangle parts are which                              
        rectcentres = [ rectcentres, [1; 0] ];
        
        spacecentres = rectcentre( [nodes(3,:); nodes(7,:)], ...
                                   [nodes(6,:); nodes(10,:)] );        
        
    elseif abs(xycoords(7,2)) <= minsizes || abs(((2*(yperiod/2)) - abs(xycoords(7,2)))) <= minsizes
        % The top of the rectangle B is against the base of the sim, or
        % equivalently the top of rectangle B is against the top of the sim
        % (or at least within 1e-4 m of it (0.1 mm)
        % Draw line 5--6 and join up nodes to top then draw bottom mag
        % rectangle 1, link nodes 4-6 and 2-bot
        
        %            8 *************** 9                                            
        %              *     8-9     *                                              
        %              *             *                                              
        %              *             *                                              
        %              *             *                                              
        %          6-8 *      B      * 7-9                                          
        %              *             *                                              
        %              *             *                                              
        %              *             *                                              
        %              *     6-7     *            
        %            6 *************** 7                                                                                          
        %              *             *                                              
        %          6-4 *             * 7-5
        %              *             *  
        %            4 *************** 5                                            
        %              *     4-5     *
        %              *             *                                              
        %              *             *               nodes:                          
        %              *             *                                              
        %          2-4 *      A      * 3-5      [ bottom left     0                                 
        %              *             *            bottom right    1                                  
        %              *             *            rectA bot left  2                                  
        %              *             *            rectA bot right 3                                  
        %              *     2-3     *            rectA top left  4                                  
        %            2 *************** 3          rectA top right 5                                  
        %              *             *            rectB bot left  6
        %          0-2 *             *  1-3       rectB bot right 7
        %              *             *            top left        8
        %            0 *************** 1          top right ]     9                                  
        %                    0-1
        
        nodes = [bottomnodes;
                 xycoords(1:4,:);
                 xycoords(5:6,:);
                 topnodes];
             
        links = [0,1; 1,3; 3,2; 2,0; 3,5; 2,4; 4,5; 5,7; 7,6; 6,4; 7,9; 6,8; 8,9];
        
        rectcentres = rectcentre( [nodes(3,:); nodes(7,:)], ...
                                  [nodes(6,:); nodes(10,:)] );    
                              
        % Add a designation to determine which rectangle parts are which                              
        rectcentres = [ rectcentres, [0; 1] ];
        
        spacecentres = rectcentre( [nodes(1,:); nodes(5,:)], ...
                                   [nodes(4,:); nodes(8,:)] );        
        
    elseif xycoords(3,2) < xycoords(1,2)
        % rectangle A is split along the top/bottom boundaries, rectangle B
        % is in middle
        % Draw line 1--2 and join up nodes to top then draw line 3--4 and
        % join up nodes to base, then draw top mag rectangle 1, link nodes
        % 4-6 and 8-2
        
        %           10 *************** 11                                           
        %              *      A      *   
        %              *             *                                                                                           
        %              *             *                                             
        %              *             *                                                                                           
        %            8 *************** 9                                                                                        
        %              *             *                                              
        %              *             *   
        %              *             *   
        %              *             *   
        %            6 *************** 7                                         
        %              *             *                                              
        %              *             *                                                                                          
        %              *             *                                              
        %              *      B      *                                                                                     
        %              *             *          nodes:                              
        %              *             *                                              
        %              *             *     [ bottom left      0
        %            4 *************** 5     bottom right     1                                                                 
        %              *             *       rectA top left   2                     
        %              *             *       rectA top right  3
        %              *             *       rectB bot left   4
        %              *             *       rectB bot right  5
        %            2 *************** 3     rectB top left   6                     
        %              *             *       rectB top right  7
        %              *             *       rectA bot left   8                                                                  
        %              *             *       rectA bot right  9                    
        %              *      A      *       top left         10                                                         
        %            0 *************** 1     top right ]      11
        %                     

        nodes = [bottomnodes;
                 xycoords(3:4,:);
                 xycoords(5:8,:);
                 xycoords(1:2,:);
                 topnodes];
             
        links = [ 0,1; 1,3; 0,2; 2,3; 3,5; 5,4; 4,2; 5,7; 4,6; 6,7; 7,9; 9,8; 8,6; 9,11; 8,10; 10,11 ];
            
        rectcentres = rectcentre( [nodes(1,:); nodes(9,:); nodes(5,:)], ...
                                  [nodes(4,:); nodes(12,:); nodes(8,:)] );
                  
        % Add a designation to determine which rectangle parts are which
        rectcentres = [rectcentres, [0; 0; 1]];
         
        spacecentres = rectcentre( [nodes(3,:); nodes(7,:)], ...
                                   [nodes(6,:); nodes(10,:)] );         
        
    elseif xycoords(7,2) < xycoords(5,2)
        % rectangle B is split along the top/bottom boundaries, rectangle A
        % is in middle
        % Draw line 7--8 and join up nodes to base then draw line 5--6 and
        % join up nodes to top, then draw bottom mag rectangle 1, link nodes
        % 4-6 and 8-2
        
        %           10 *************** 11                                           
        %              *      B      *   
        %              *             *                                                                                           
        %              *             *                                             
        %              *             *                                                                                           
        %            8 *************** 9                                                                                        
        %              *             *                                              
        %              *             *   
        %              *             *   
        %              *             *   
        %            6 *************** 7                                         
        %              *             *                                              
        %              *             *                                                                                          
        %              *             *                                              
        %              *      A      *                                                                                     
        %              *             *          nodes:                              
        %              *             *                                              
        %              *             *     [ bottom left      0
        %            4 *************** 5     bottom right     1                                                                 
        %              *             *       rectB top left   2                     
        %              *             *       rectB top right  3
        %              *             *       rectA bot left   4
        %              *             *       rectA bot right  5
        %            2 *************** 3     rectA top left   6                     
        %              *             *       rectA top right  7
        %              *             *       rectB bot left   8                                                                  
        %              *             *       rectB bot right  9                    
        %              *      B      *       top left         10                                                         
        %            0 *************** 1     top right ]      11
        %                     

        nodes = [bottomnodes;
                 xycoords(7:8,:);
                 xycoords(1:4,:);
                 xycoords(5:6,:);
                 topnodes];
             
        links = [ 0,1; 1,3; 0,2; 2,3; 3,5; 5,4; 4,2; 5,7; 4,6; 6,7; 7,9; 9,8; 8,6; 9,11; 8,10; 10,11 ];        
        
        rectcentres = rectcentre( [nodes(1,:); nodes(9,:); nodes(5,:)], ...
                                  [nodes(4,:); nodes(12,:); nodes(8,:)] );

        % Add a designation to determine which rectangle parts are which                              
        rectcentres = [rectcentres, [1; 1; 0]];
         
        spacecentres = rectcentre( [nodes(3,:); nodes(7,:)], ...
                                   [nodes(6,:); nodes(10,:)] );         
         
    elseif xycoords(1,2) < xycoords(3,2) && xycoords(5,2) < xycoords(7,2) && xycoords(5,2) > xycoords(1,2)
        % Draw both rectangles, rectangle A at bottom and rectangle B at
        % top then link nodes 4-6 and 8-top and 2-bot
        
        %
        %           10 *************** 11                                         
        %              *             *                                              
        %              *             * 
        %            8 *************** 9                                          
        %              *             *                                              
        %              *             *                                                                                          
        %              *             *                                              
        %              *      B      *                                                                                     
        %              *             *                                              
        %              *             *                                              
        %              *             *            
        %            6 *************** 7         nodes:                                                                         
        %              *             *                                              
        %              *             *       [ bottom left     0
        %              *             *         bottom right    1 
        %              *             *         rectA bot left  2
        %            4 *************** 5       rectA bot right 3                    
        %              *             *         rectA top left  4
        %              *             *         rectA top right 5                                                                 
        %              *             *         rectB bot left  6                   
        %              *      A      *         rectB bot left  7                               
        %              *             *         rectB top left  8                                     
        %              *             *         rectB top left  9
        %              *             *         top left        10                 
        %            2 *************** 3       top right ]     11                  
        %              *             *  
        %              *             *  
        %            0 *************** 1                                          
        %              
        %   
        
        nodes = [ bottomnodes;
                  xycoords(1:8,:);
                  topnodes ];
              
        links = [ 0,1; 1,3; 0,2; 2,3; 3,5; 5,4; 4,2; 5,7; 4,6; 6,7; 7,9; 9,8; 8,6; 9,11; 8,10; 10,11 ];
        
        rectcentres = rectcentre( [nodes(3,:); nodes(7,:)], ...
                                  [nodes(6,:); nodes(10,:)] );    
                   
        % Add a designation to determine which rectangle parts are which
        rectcentres = [ rectcentres, [0; 1] ];
        
        spacecentres = rectcentre( [nodes(1,:); nodes(5,:); nodes(9,:)], ...
                                   [nodes(4,:); nodes(8,:); nodes(12,:)] );        
        
    elseif xycoords(1,2) < xycoords(3,2) && xycoords(5,2) < xycoords(7,2) && xycoords(5,2) < xycoords(1,2)
        % Draw both rectangles, rectangle A at top and rectangle B at bottom
        % then link nodes 2-8 and 4-top and 6-bot
        
        %
        %           10 *************** 11                                         
        %              *             *                                              
        %              *             * 
        %            8 *************** 9                                          
        %              *             *                                              
        %              *             *                                                                                          
        %              *             *                                              
        %              *      A      *                                                                                     
        %              *             *                                              
        %              *             *                                              
        %              *             *            
        %            6 *************** 7         nodes:                                                                         
        %              *             *                                              
        %              *             *       [ bottom left     0
        %              *             *         bottom right    1 
        %              *             *         rectB bot left  2
        %            4 *************** 5       rectB bot right 3                    
        %              *             *         rectB top left  4
        %              *             *         rectB top right 5                                                                 
        %              *             *         rectA bot left  6                   
        %              *      B      *         rectA bot left  7                               
        %              *             *         rectA top left  8                                     
        %              *             *         rectA top left  9
        %              *             *         top left        10                 
        %            2 *************** 3       top right ]     11                  
        %              *             *  
        %              *             *  
        %            0 *************** 1                                          
        %              
        %   
        
        nodes = [ bottomnodes;
                  xycoords(5:8,:);
                  xycoords(1:4,:);
                  topnodes ];
              
        links = [ 0,1; 1,3; 0,2; 2,3; 3,5; 5,4; 4,2; 5,7; 4,6; 6,7; 7,9; 9,8; 8,6; 9,11; 8,10; 10,11 ];
        
        rectcentres = rectcentre( [nodes(3,:); nodes(7,:)], ...
                                  [nodes(6,:); nodes(10,:)] );    
                              
        % Add a designation to determine which rectangle parts are which
        rectcentres = [ rectcentres, [1; 0] ];
        
        spacecentres = rectcentre( [nodes(1,:); nodes(5,:); nodes(9,:)], ...
                                   [nodes(4,:); nodes(8,:); nodes(12,:)] );           
        
    end
    
    % add the x offset to the x-coords
    nodes(:,1) = nodes(:,1) + xoffset;
    rectcentres(:,1) = rectcentres(:,1) + xoffset;
    spacecentres(:,1) = spacecentres(:,1) + xoffset;
    
    % add a count to the nodes
    nodeids = (0:size(nodes, 1)-1) + nodecount;
    
    % add the starting link number to the links
    links = links + nodecount;

end