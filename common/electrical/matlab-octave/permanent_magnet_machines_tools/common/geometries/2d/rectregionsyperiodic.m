function [nodes, nodeids, links, rectcentres, spacecentres] = ...
    rectregionsyperiodic(x, y1, y2, xoffset, ydisp, varargin)
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
% rectregionsyperiodic also accepts several optionsl arguments as
% parameter-value pairs:
%
%  'Tol' - optional tolerance at which distances bewteen elements are
%    considered to be of size zero. If not supplied, 1e-5 is used.
%
%  'NodeCount' - optional starting count for the node ids used to specify the
%    links between nodes. This number will be added to all the link
%    numbers. Defaults to zero if not supplied.  
%
%  'NY1Pairs' - optional number of pairs of regions to create, default is 1
%
% Output
%
%  nodes - (n x 2) matrix of x and y coordinates of the geometry. The nodes are
%    added in sequence starting from the bottom left node, followed by the 
%    bottom right node, then the next lowest left node and right node and so on
%
%  nodeids -  vector of id numbers for the nodes in nodes
%
%  links - (m x 2) matrix of links bewteen nodes specified as the two id
%    numbers of the nodes. The first link is always the horizontal link between
%    nodes 1 and 2 at the bottom of the geometry, and the last link is always
%    the horizontal link between the final two nodes at the top of the geometry.
%
%  rectcentres - (p x 2) matrix of x and y coordinates of the centers of
%    the 'A' and 'B' rectangles, or their component parts (can be split
%    into two when 'wrapped')
%
%  spacecentres - (p x 2) matrix of x and y coordinates of the centers of
%    the spaces between the 'A' and 'B' rectangles.
%
%

    Inputs.Tol = 1e-5;
    Inputs.NodeCount = 0;
    Inputs.NY1Pairs = 1;
    
    Inputs = parse_pv_pairs (Inputs, varargin);
    
    Ny1 = 2 * Inputs.NY1Pairs;
    
    yperiod = 2*y1 + 2*y2;
    
    % First define the rectangle 1 coordinates, the other rectangles will
    % be defined relative to these
    % Bottom left of bottom rect
    xycoords(1,:) = [-x/2, y2/2];
    % Bottom right of bottom rect
    xycoords(2,:) = [x/2, y2/2];
    % Top left of bottom rect
    xycoords(3,:) = [-x/2, xycoords(1,2) + y1];
    % Top right of bottom rect
    xycoords(4,:) = [x/2, xycoords(2,2) + y1];
    % Next rectangles are bottom rectangle shifted by n half periods
    
    for n = 1:(Ny1-1)
        xycoords = [ xycoords; 
                     xycoords(1:4,1), xycoords(end-3:end,2) + yperiod/2; ];
    end
    
    % Label Coordinates, bottom rect then top rect
    xycoords = [ xycoords;
                 (xycoords(1:4:end,1)+xycoords(2:4:end,1))/2, (xycoords(1:4:end,2)+xycoords(3:4:end,2))/2; ];
    
%     xycoords(9,:) = [xycoords(1,1)+(xycoords(2,1)-xycoords(1,1))/2, xycoords(1,2)+(xycoords(3,2)-xycoords(1,2))/2];
%     xycoords(10,:) = [xycoords(5,1)+(xycoords(6,1)-xycoords(5,1))/2, xycoords(5,2)+(xycoords(7,2)-xycoords(5,2))/2];
    
    % Map coordinates onto cylinder surface
    wrapperiod = (yperiod/2) * Ny1;
    xycoords(:,2) = xycoords(:,2) .* pi ./ wrapperiod;
    
    wrapydisp = ydisp .* pi ./ wrapperiod;
    
    % rotate around cylinder
    xycoords(:,2) = rem(xycoords(:,2) + wrapydisp, pi);
    
    % negative coordinates imply they are in fact rotated backwards which
    % we will take to be positions relative to the top of the sim
    xycoords(sign(xycoords(:,2)) == -1, 2) = pi + xycoords(sign(xycoords(:,2)) == -1, 2); 
     
    % Map back to x-y plane
    xycoords(:,2) = wrapperiod .* xycoords(:,2) ./ pi;
    
    minsizes = max(y1 * 0.0005, Inputs.Tol);
    
    bottomnodes = [-x/2, 0; x/2, 0];
    topnodes = [-x/2, wrapperiod; x/2, wrapperiod];
    
    % now generate the nodes and joints according to the specified
    % tolerances
    if any( abs(xycoords(1:8:Ny1*4,2)) <= minsizes ...
            | abs(((Ny1*(yperiod/2)) - abs(xycoords(1:8:Ny1*4,2)))) <= minsizes )
        
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
        
        Apos = find ( abs(xycoords(1:8:Ny1*4,2)) <= minsizes ...
                       | abs(((Ny1*(yperiod/2)) - abs(xycoords(1:8:Ny1*4,2)))) <= minsizes ) * 2 - 1;
                   
        tempnodes = circshift (xycoords(1:4*Ny1,:), [-Apos*4+2, 0]);
        
        nodes = [bottomnodes; 
                 tempnodes(1:end-2,:);
                 ... xycoords(3:4*Ny1,:);
                 ... xycoords(5:8,:); 
                 topnodes];
        
        links = [ (0:2:4*Ny1)', (0:2:4*Ny1)'+1;
                  (0:(4*Ny1-1))', (0:(4*Ny1-1))'+2 ];
%         links = [0,1; 2,3; 4,5; 6,7; 8,9; 0,2; 1,3; 2,4; 3,5; 4,6; 5,7; 6,8; 7,9;]; 

%         rectcentres = rectcentre( [nodes(1,:); nodes(5,:)], ...
%                                   [nodes(4,:); nodes(8,:)] );
        rectcentres = [ rectcentre( nodes(1:4:((4*Ny1)-3),:), ...
                                    nodes(4:4:4*Ny1,:) ), ...
                        nan * ones(Ny1, 1) ];
                              
        % Add a designation to determine which rectangle parts are which
        desig = [0,1];
        for n = 1:Ny1
            rectcentres(n,3) = desig(1);
            desig = fliplr (desig);
        end
        
        spacecentres = rectcentre( nodes(3:4:((4*Ny1)-1),:), ...
                                   nodes(6:4:4*Ny1+2,:) );
%         spacecentres = rectcentre( [nodes(3,:); nodes(7,:)], ...
%                                    [nodes(6,:); nodes(10,:)] );

    elseif any( abs(xycoords(3:8:4*Ny1+2,2)) <= minsizes ...
                | abs(((Ny1*(yperiod/2)) - abs(xycoords(3:8:4*Ny1+2,2)))) <= minsizes )
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

        Apos = find ( abs(xycoords(3:8:4*Ny1+2,2)) <= minsizes ...
                       | abs(((Ny1*(yperiod/2)) - abs(xycoords(3:8:4*Ny1+2,2)))) <= minsizes ) * 2 - 1;
                   
        tempnodes = circshift (xycoords(1:4*Ny1,:), [-Apos*4, 0]);
        
        nodes = [bottomnodes; 
                 tempnodes(1:end-2,:);
                 ... xycoords(5:4*Ny1,:); 
                 ... xycoords(1:2,:); 
                 topnodes];
             
        links = [ (0:2:4*Ny1)', (0:2:4*Ny1)'+1;
                  (0:(4*Ny1-1))', (0:(4*Ny1-1))'+2 ];
%         links = [0,1; 2,3; 0,2; 1,3; 3,5; 5,4; 4,2; 4,6; 5,7; 6,7; 6,8; 7,9; 8,9]; 

%         rectcentres = rectcentre( [nodes(3,:); nodes(7,:)], ...
%                                   [nodes(6,:); nodes(10,:)] );
        rectcentres = [ rectcentre( nodes(3:4:((4*Ny1)-1),:), ...
                                    nodes(6:4:4*Ny1+2,:) ), ...
                        nan * ones(Ny1, 1) ];
                  
        % Add a designation to determine which rectangle parts are which
%         rectcentres = [ rectcentres, [1; 0] ];
        desig = [1,0];
        for n = 1:Ny1
            rectcentres(n,3) = desig(1);
            desig = fliplr (desig);
        end
        
%         spacecentres = rectcentre( [nodes(1,:); nodes(5,:)], ...
%                                    [nodes(4,:); nodes(8,:)] );
        spacecentres = rectcentre( nodes(1:4:((4*Ny1)-3),:), ...
                                   nodes(4:4:4*Ny1,:) );
        
    elseif any( abs(xycoords(5:8:4*Ny1-3,2)) <= minsizes ...
                | abs(((Ny1*(yperiod/2)) - abs(xycoords(5:8:4*Ny1-3,2)))) <= minsizes )
        
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

        Bpos = find ( abs(xycoords(5:8:4*Ny1-3,2)) <= minsizes ...
                      | abs(((Ny1*(yperiod/2)) - abs(xycoords(5:8:4*Ny1-3,2)))) <= minsizes ...
                    ) * 2;
        
        tempnodes = circshift (xycoords(1:4*Ny1,:), [-Bpos*4+2, 0]);
        
%         nodes = [bottomnodes; 
%                  xycoords(7:8,:); 
%                  xycoords(1:4,:); 
%                  topnodes];
        nodes = [bottomnodes; 
                 tempnodes(1:end-2,:);
                 ...xycoords(((Ny1*4)-1):(Ny1*4),:); 
                 ...xycoords(1:((Ny1-1)*4),:); 
                 topnodes];
             
        links = [ (0:2:4*Ny1)', (0:2:4*Ny1)'+1;
                  (0:(4*Ny1-1))', (0:(4*Ny1-1))'+2 ];
%         links = [0,1; 2,3; 4,5; 6,7; 8,9; 0,2; 1,3; 2,4; 3,5; 4,6; 5,7; 6,8; 7,9];
        
%         rectcentres = rectcentre( [nodes(1,:); nodes(5,:)], ...
%                                   [nodes(4,:); nodes(8,:)] );
        rectcentres = [ rectcentre( nodes(1:4:((4*Ny1)-3),:), ...
                                    nodes(4:4:4*Ny1,:) ), ...
                        nan * ones(Ny1, 1) ];
                              
        % Add a designation to determine which rectangle parts are which
        desig = [1,0];
        for n = 1:Ny1
            rectcentres(n,3) = desig(1);
            desig = fliplr (desig);
        end
              
%         % Add a designation to determine which rectangle parts are which                              
%         rectcentres = [ rectcentres, [1; 0] ];
        
%         spacecentres = rectcentre( [nodes(3,:); nodes(7,:)], ...
%                                    [nodes(6,:); nodes(10,:)] );        
        spacecentres = rectcentre( nodes(3:4:((4*Ny1)-1),:), ...
                                   nodes(6:4:4*Ny1+2,:) );
        
    elseif any(abs(xycoords(7:8:4*Ny1-1,2)) <= minsizes ...
            | abs(((Ny1*(yperiod/2)) - abs(xycoords(7:8:4*Ny1-1,2)))) <= minsizes )
        
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
        
        
        Bpos = find ( abs(xycoords(7:8:4*Ny1-1,2)) <= minsizes ...
                      | abs(((Ny1*(yperiod/2)) - abs(xycoords(7:8:4*Ny1-1,2)))) <= minsizes ...
                    ) * 2;
        
        tempnodes = circshift (xycoords(1:4*Ny1,:), [-Bpos*4, 0]);
        
        nodes = [bottomnodes;
                 tempnodes(1:end-2,:);
                 ... xycoords(5:6,:);
                 topnodes];
             
%         links = [0,1; 2,3; 4,5; 6,7; 8,9;  0,2; 1,3; 2,4; 3,5; 4,6; 5,7; 6,8; 7,9;];
        links = [ (0:2:4*Ny1)', (0:2:4*Ny1)'+1;
                  (0:(4*Ny1-1))', (0:(4*Ny1-1))'+2 ];
        
%         rectcentres = rectcentre( [nodes(3,:); nodes(7,:)], ...
%                                   [nodes(6,:); nodes(10,:)] );    
%                               
%         % Add a designation to determine which rectangle parts are which                              
%         rectcentres = [ rectcentres, [0; 1] ];
%         
%         spacecentres = rectcentre( [nodes(1,:); nodes(5,:)], ...
%                                    [nodes(4,:); nodes(8,:)] );
                               
        rectcentres = [ rectcentre( nodes(3:4:((4*Ny1)-1),:), ...
                                    nodes(6:4:4*Ny1+2,:) ), ...
                        nan * ones(Ny1, 1) ];
                  
        % Add a designation to determine which rectangle parts are which
%         rectcentres = [ rectcentres, [1; 0] ];
        desig = [0,1];
        for n = 1:Ny1
            rectcentres(n,3) = desig(1);
            desig = fliplr (desig);
        end
        
%         spacecentres = rectcentre( [nodes(1,:); nodes(5,:)], ...
%                                    [nodes(4,:); nodes(8,:)] );
        spacecentres = rectcentre( nodes(1:4:((4*Ny1)-3),:), ...
                                   nodes(4:4:4*Ny1,:) );
        
    elseif any(xycoords(3:8:4*Ny1+2,2) < xycoords(1:8:4*Ny1,2))
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

        Apos = find (xycoords(3:8:4*Ny1+2,2) < xycoords(1:8:4*Ny1,2)) * 2 - 1;
        
        tempnodes = circshift (xycoords(1:4*Ny1,:), [-Apos*4+2, 0]);
        
        nodes = [bottomnodes;
                 tempnodes;
                 ...xycoords(Apos*4-1:4*Ny1,:); 
                 ... xycoords(3:4,:);
                 ... xycoords(5:8,:);
                 ...xycoords(1:Apos*4-2,:);
                 topnodes];
             
%         links = [ 0,1; 2,3; 4,5; 6,7; 8,9; 10,11; 0,2; 1,3; 2,4; 3,5; 4,6; 5,7; 6,8; 7,9; 8,10; 9,11; ];
        links = [ (0:2:4*(Ny1+1)-1)', (0:2:4*(Ny1+1)-1)'+1;
                  (0:(4*Ny1+1))', (0:(4*Ny1+1))'+2 ];
              
        rectcentres = [ rectcentre( [nodes(1,:); nodes((Ny1*4)+1,:); nodes(5:4:((Ny1-1)*4+1),:)], ...
                                    [nodes(4,:); nodes((Ny1*4)+4,:); nodes(8:4:((Ny1-1)*4+4),:)] ), ...
                        nan * ones(Ny1+1, 1) ];
                  
        % Add a designation to determine which rectangle parts are which
%         rectcentres = [rectcentres, [0; 0; 1]];
        rectcentres(1,3) = 0;
        rectcentres(2,3) = 0;
        desig = [1,0];
        for n = 3:Ny1+1
            rectcentres(n,3) = desig(1);
            desig = fliplr (desig);
        end
         
%         spacecentres = rectcentre( [nodes(3,:); nodes(7,:)], ...
%                                    [nodes(6,:); nodes(10,:)] );
        spacecentres = rectcentre( nodes(3:4:((4*Ny1)-1),:), ...
                                   nodes(6:4:4*Ny1+2,:) );
        
    elseif any(xycoords(7:8:4*Ny1-1,2) < xycoords(5:8:4*Ny1-3,2))
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

        Bpos = find ( xycoords(7:8:4*Ny1-1,2) < xycoords(5:8:4*Ny1-3,2) ) * 2;
        
        tempnodes = circshift (xycoords(1:4*Ny1,:), [-Bpos*4+2, 0]);
        
        nodes = [bottomnodes;
                 tempnodes;
                 ...xycoords(Ny1*4-1:Ny1*4,:);
                 ...xycoords(1:Ny1*4-2,:);
                 ... xycoords(7:8,:);
                 ... xycoords(1:4,:);
                 ... xycoords(5:6,:);
                 topnodes];
             
%         links = [ 0,1; 1,3; 0,2; 2,3; 3,5; 5,4; 4,2; 5,7; 4,6; 6,7; 7,9; 9,8; 8,6; 8,10; 9,11; 10,11 ];        
%         
%         rectcentres = rectcentre( [nodes(1,:); nodes(9,:); nodes(5,:)], ...
%                                   [nodes(4,:); nodes(12,:); nodes(8,:)] );
                              
        links = [ (0:2:4*(Ny1+1)-1)', (0:2:4*(Ny1+1)-1)'+1;
                  (0:(4*Ny1+1))', (0:(4*Ny1+1))'+2 ];
              
        rectcentres = [ rectcentre( [nodes(1,:); nodes((Ny1*4)+1,:); nodes(5:4:((Ny1-1)*4+1),:)], ...
                                    [nodes(4,:); nodes((Ny1*4)+4,:); nodes(8:4:((Ny1-1)*4+4),:)] ), ...
                        nan * ones(Ny1+1, 1) ];
                    
        % Add a designation to determine which rectangle parts are which                              
%         rectcentres = [rectcentres, [1; 1; 0]];
        rectcentres(1,3) = 1;
        rectcentres(2,3) = 1;
        desig = [0,1];
        for n = 3:Ny1+1
            rectcentres(n,3) = desig(1);
            desig = fliplr (desig);
        end
         
%         spacecentres = rectcentre( [nodes(3,:); nodes(7,:)], ...
%                                    [nodes(6,:); nodes(10,:)] );         
        spacecentres = rectcentre( nodes(3:4:((4*Ny1)-1),:), ...
                                   nodes(6:4:4*Ny1+2,:) );
         
    else
        
        remydisp = abs(rem (ydisp, yperiod));
% unwind even number of turns and then look at node locations
        if ( ydisp == 0 ) ...
                || (ydisp > 0 && ( remydisp < (yperiod/4) || remydisp >= (3*yperiod/4)) ) ...
                || (ydisp < 0 && ( remydisp <= (yperiod/4) || remydisp > (3*yperiod/4)) )
        
            % First A is below first B
            %
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

            [~,Apos] = min (xycoords(1:4*Ny1,2));

            tempnodes = circshift (xycoords(1:4*Ny1,:), [-(Apos-1), 0]);

            nodes = [ bottomnodes;
                      tempnodes;
                      ...xycoords(1:4*Ny1,:);
                      ... xycoords(1:8,:);
                      topnodes ];

    %         links = [ 0,1; 1,3; 0,2; 2,3; 3,5; 5,4; 4,2; 5,7; 4,6; 6,7; 7,9; 9,8; 8,6; 8,10; 9,11; 10,11 ];
            links = [ (0:2:4*(Ny1+1)-1)', (0:2:4*(Ny1+1)-1)'+1;
                      (0:(4*Ny1+1))', (0:(4*Ny1+1))'+2 ];

    %         rectcentres = rectcentre( [nodes(3,:); nodes(7,:)], ...
    %                                   [nodes(6,:); nodes(10,:)] );    
    %                    
    %         % Add a designation to determine which rectangle parts are which
    %         rectcentres = [ rectcentres, [0; 1] ];
            rectcentres = [ rectcentre( nodes(3:4:((4*Ny1)-1),:), ...
                                        nodes(6:4:4*Ny1+2,:) ), ...
                            nan * ones(Ny1, 1) ];

            % Add a designation to determine which rectangle parts are which
    %         rectcentres = [ rectcentres, [1; 0] ];
            desig = [0,1];
            for n = 1:Ny1
                rectcentres(n,3) = desig(1);
                desig = fliplr (desig);
            end

    %         spacecentres = rectcentre( [nodes(1,:); nodes(5,:); nodes(9,:)], ...
    %                                    [nodes(4,:); nodes(8,:); nodes(12,:)] );   
            spacecentres = rectcentre( nodes(1:4:4*Ny1+2,:), nodes(4:4:4*(Ny1+1),:) );

        else
            % First A is above first B
            %
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

            [~,Apos] = min (xycoords(1:4*Ny1,2));

            tempnodes = circshift (xycoords(1:4*Ny1,:), [-(Apos-1), 0]);

            nodes = [ bottomnodes;
                      tempnodes;
                      ...xycoords(5:4*Ny1,:);
                      ... xycoords(5:8,:);
                      ... xycoords(1:4,:);
                      topnodes ];

    %         links = [ 0,1; 1,3; 0,2; 2,3; 3,5; 5,4; 4,2; 5,7; 4,6; 6,7; 7,9; 9,8; 8,6; 8,10; 9,11; 10,11 ];
            links = [ (0:2:4*(Ny1+1)-1)', (0:2:4*(Ny1+1)-1)'+1;
                      (0:(4*Ny1+1))', (0:(4*Ny1+1))'+2 ];
    %         rectcentres = rectcentre( [nodes(3,:); nodes(7,:)], ...
    %                                   [nodes(6,:); nodes(10,:)] );    
    %                               
    %         % Add a designation to determine which rectangle parts are which
    %         rectcentres = [ rectcentres, [1; 0] ];
            rectcentres = [ rectcentre( nodes(3:4:((4*Ny1)-1),:), ...
                                        nodes(6:4:4*Ny1+2,:) ), ...
                            nan * ones(Ny1, 1) ];

            desig = [1,0];
            for n = 1:Ny1
                rectcentres(n,3) = desig(1);
                desig = fliplr (desig);
            end

%             spacecentres = rectcentre( [nodes(1,:); nodes(5,:); nodes(9,:)], ...
%                                        [nodes(4,:); nodes(8,:); nodes(12,:)] );           

            spacecentres = rectcentre( nodes(1:4:4*Ny1+2,:), nodes(4:4:4*(Ny1+1),:) );           

        end
    end
    
    % add the x offset to the x-coords
    nodes(:,1) = nodes(:,1) + xoffset;
    rectcentres(:,1) = rectcentres(:,1) + xoffset;
    spacecentres(:,1) = spacecentres(:,1) + xoffset;
    
    % add a count to the nodes
    nodeids = (0:size(nodes, 1)-1) + Inputs.NodeCount;
    
    % add the starting link number to the links
    links = links + Inputs.NodeCount;

end