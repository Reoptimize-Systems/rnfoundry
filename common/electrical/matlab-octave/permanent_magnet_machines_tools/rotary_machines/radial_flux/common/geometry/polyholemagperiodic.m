function [FemmProblem, nodes, nodeids, links, magblockinds] = ...
    polyholemagperiodic(FemmProblem, thetapole, thetamag, rmag, roffset, pos, varargin)
% generates a periodic magnets containing region with straight sides forming a 
% closed semi-regular polygon shape (i.e. a polygon with exactly two side 
% lengths, with a hole the smae shape) suitible for modelling radial flux
% machines with rectangular (rather than curved) magnets.
%
% Syntax
%
% [FemmProblem, nodes, nodeids, links] = ...
%           polyholemagperiodic(FemmProblem, thetapole, thetamag, rmag, roffset, pos)
% [...] = polyholemagperiodic(..., 'Paramter', Value)
%
%
% Description
%
% polyholemagperiodic creates a periodic geometry of two magnets with
% spaces in between, with a base position shown in the figure below.
%
%             ********
%       ******        *                                    
%         *             *                                              
%           *        ****** 
%            *********     *                                            
%             *             *                                              
%               *             *                                                                                          
%                *             *                                              
%                  *    Mag 2    *                                                                                     
%                   *             *                                              
%                    *             *                                              
%                     *             *            
%                      ***************                                                                                          
%                       *             *                                              
%                        *             *  .........................
%                         *             *                ^
%                         *             *                 :
%                          *************** ..^.......      :                            
%                          *             *   :             :
%                           *             *   :             :                                                                          
%                           *             *   : thetamag    :                                    
%                           *    Mag 1    *    :             : thetapole             
%                            *             *   :             :                                      
%                            *             *   :             :            
%                            *             * . v.........     :                                       
%                             ***************                 :                                       
%                             *             *                 :               
%                             *             *                 v                
%  x                          *************** ..............................                                       
% r=0                         <------------->
%  :                               rmag
%  :                                :
%  :                                :
%  :            roffset             :
%  :------------------------------->:
%  :           
%
%
% This geometry is drawn in a periodic way in the tangential direction,
% 'wrapping' around at the top and bottom. 
%
% Inputs
%
%  FemmProblem - FemmProblem structure to which the geometry will be added
%
%  thetapole - pole width in radians
%
%  thetamag - magnet width in radians
%
%  rmag - radial thickness of the magnets
%
%  roffset - radial displacement of the magnet centers from the center
%
%  pos - the angular position of the magnets
%
%  In addition, a number of optional parameters can be specified as
%  parameter-value pairs. Possible parameter-value pairs are:
%
%  'MagDirections' - either a 2 element numeric vector, or a 2 element cell
%    array of strings. If numeric, these are the directions in degrees of
%    the magnet magnetisation. If a cell array of strings, these are
%    evaluated in the FEMM or xfemm lua interpreter to yield the magnet
%    direction in the magnet region elements. Variables that can be used in
%    these strings are:
%
%    'theta': angle in degrees of a line connecting the center of each
%             element with the origin 
%
%    'R'    : length of a line connecting the center of each element with the
%             origin
%
%    'x'    : x position of each element
%
%    'y'    : y position of each elements
% 
%    The default is {'theta', 'theta+180'}, resulting in radially
%    magnetized magnets of opposite polarity.
%
%  'MagnetMaterial' - 
%
%  'MagnetGroup' - 
%
%  'SpaceMaterial' - 
%
%  'SpaceGroup' - 
%
%  'Tol' - 
%
%  'MeshSize' - 
%
% Output
%
%  FemmProblem - 
%
%  nodes - 
%
%  nodeids - 
%
%  links - 
%
%  magblockinds - 
%
%

    Inputs.MagDirections = {'theta', 'theta+180'};
    Inputs.MagnetMaterial = 1;
    Inputs.MagnetGroup = 0;
    Inputs.SpaceMaterial = 1;
    Inputs.SpaceGroup = 0;
    Inputs.Tol = 1e-5;
    Inputs.MeshSize = -1;
    Inputs.NPolePairs = 1;

    Inputs = parse_pv_pairs(Inputs, varargin);
    
    if isnumeric(Inputs.MagDirections) && isvector(Inputs.MagDirections)
        Inputs.MagDirections = {Inputs.MagDirections(1), Inputs.MagDirections(2)};
    end
    
    % get the number of existing nodes, segments, boundaries etc. if any
    elcount = elementcount_mfemm(FemmProblem);
    
    % add a periodic boundary for the edges of the magnets region
    [FemmProblem, boundind] = addboundaryprop_mfemm(FemmProblem, ...
                                                    'Rect Mags Periodic', 4);
    
    % construct the segments and label locations for the magnet regions
    % first in a linear fashion which we will manipulate into the real
    % shape by modifying the node locations
    [nodes, nodeids, links, rectcentres, spacecentres] = ...
        rectregionsyperiodic(rmag, thetamag, (thetapole-thetamag), roffset, pos, ...
                             'Tol', Inputs.Tol, ...
                             'NodeCount', elcount.NNodes, ...
                             'NY1Pairs', Inputs.NPolePairs);

    % get the vertical links by finding those links where the difference in
    % y coordinates of the link nodes is not zero, these links must be made
    % into arc segments
%    vertlinks = links(abs(diff( [nodes(links(:,1)+1-elcount.NNodes,1), nodes(links(:,2)+1-elcount.NNodes,1)], 1, 2 )) < Inputs.Tol,:);
%    angles = diff( [nodes(vertlinks(:,1)+1-elcount.NNodes,2), nodes(vertlinks(:,2)+1-elcount.NNodes,2)], 1, 2);
    
    % get the horizontal links, these will be segments
%    horizlinks = links(abs(diff( [nodes(links(:,1)+1-elcount.NNodes,1), nodes(links(:,2)+1-elcount.NNodes,1)], 1, 2 )) >= Inputs.Tol,:);
    
    xmid = mean([nodes(1,1), nodes(2,1)]);
    
    allmidpoints = [ repmat(xmid, size(nodes,1)/2, 1), nodes(1:2:end-1,2) ];
    
    % transform the node locations to put the nodes onto the circumference of a
    % circle
    [allmidpoints(:,1), allmidpoints(:,2)] = pol2cart(allmidpoints(:,2), allmidpoints(:,1));
    
    % get vectors pointing from midpoint to midpoint
    allmidpntvecs = allmidpoints(2:end,:) - allmidpoints(1:end-1,:);
    
    % 
    [transnodes(:,1), transnodes(:,2)] = pol2cart(nodes(:,2), nodes(:,1));
    
    rectcentres = sortrows (rectcentres, 2);
    spacecentres = sortrows (spacecentres, 2);
    
    thetaspace = thetapole-thetamag;
    
    % flag for indicating if the top and bottom edges of the geometry have been
    % modified
    edgemods = false;
    
    % create a new set of nodes which make rectangles with sides parallel to the
    % midpoint vectors in magnet rectangles
    for i = 1:size (rectcentres,1)
    
      segaboveind = find (nodes(:,2) > rectcentres(i,2));
      segbelowind = find (nodes(:,2) < rectcentres(i,2));
      
      % the first node found above the center will be the left hand node
      segaboveind = segaboveind(1);
      % the last node found below the center will be the right hand node, so get
      % the ind before to get the left hand node
      segbelowind = segbelowind(end-1);
      
      if size (rectcentres,1) == 2
          % the magnets are not split in the current position, the geometry is
          % split at the spacers
      
          % get mid points of horizontal segments above and below
          midpoints = [ xmid, nodes(segaboveind,2); ...
                        xmid, nodes(segbelowind,2); ];
          
          % transform to lie on arc
          [midpoints(:,1), midpoints(:,2)] = pol2cart(midpoints(:,2), midpoints(:,1));
          
          % get vector between midpoints pointing from lower to upper
          midpntvec = midpoints(1,:) - midpoints(2,:);
          
          % get new node locations of magnets on left and right
          leftnodes = bsxfun(@plus, midpoints, rmag/2 * unit ([-midpntvec(:,2), midpntvec(:,1)]));
          rightnodes = bsxfun(@plus, midpoints, rmag/2 * unit ([midpntvec(:,2), -midpntvec(:,1)]));
          
          if i == 1
          
              if rectcentres(1,2) > spacecentres(1,2)
          
                  % magnet not resting on bottom
                  spacesegaboveind = find (nodes(:,2) > spacecentres(i,2));
                  spacesegaboveind = spacesegaboveind(1);
                  
                  % get mid points of the horizontal segments above spacer label
                  midpoints = [ xmid, nodes(spacesegaboveind,2); ];
                  
                  % transform to lie on arc
                  [midpoints(:,1), midpoints(:,2)] = pol2cart(midpoints(:,2), midpoints(:,1));
                  
                  % get the midpoint of the other edge of the magnet by rotating
                  % the top mipoint about (0,0) by thetamag radians (anticlockwise
                  % with this formulation)
                  rotM = [ cos(thetaspace)  -sin(thetaspace); 
                           sin(thetaspace)   cos(thetaspace) ];
                           
                  midpoints = [ midpoints; 
                                midpoints * rotM ];
                                
                  % get vector between midpoints pointing from lower to upper
                  midpntvec = midpoints(1,:) - midpoints(2,:);
                  
                  % shift the lower midpoint up to lie at the split point of the magnet
                  midpoints(2,:) = midpoints (1,:) - midpntvec * (nodes(segaboveind,2) / thetaspace);
                  
                  % for the bottom segment nodes, at the boundary of the region, we
                  % need to find the intercept of the sides of the geometry and a radial
                  % line going from the origin
                  xy1 = [ transnodes(spacesegaboveind,:), transnodes(spacesegaboveind,:) - 10*midpntvec; ...
                          transnodes(spacesegaboveind+1,:), transnodes(spacesegaboveind+1,:) - 10*midpntvec; ...
                          0, 0, 0, 0;
                          leftnodes(2,:), rightnodes(2,:) ];
                          
              else
                  % magnet resting on bottom
                  
                  % get mid points of the horizontal segments above 
                  midpoints = [ xmid, nodes(segaboveind,2); ];
                  
                  % transform to lie on arc
                  [midpoints(:,1), midpoints(:,2)] = pol2cart(midpoints(:,2), midpoints(:,1));
                  
                  % get the midpoint of the other edge of the magnet by rotating
                  % the top mipoint about (0,0) by thetamag radians (anticlockwise
                  % with this formulation)
                  rotM = [ cos(thetamag)  -sin(thetamag); 
                           sin(thetamag)   cos(thetamag) ];
                           
                  midpoints = [ midpoints; 
                                midpoints * rotM ];
                                
                  % get vector between midpoints pointing from lower to upper
                  midpntvec = midpoints(1,:) - midpoints(2,:);
                  
                  % shift the lower midpoint up to lie at the split point of the magnet
                  midpoints(2,:) = midpoints (1,:) - midpntvec * (nodes(segaboveind,2) / thetamag);
                  
                  % get new node locations of magnets on left and right
                  leftnodes = bsxfun(@plus, midpoints, rmag/2 * unit ([-midpntvec(:,2), midpntvec(:,1)]));
                  rightnodes = bsxfun(@plus, midpoints, rmag/2 * unit ([midpntvec(:,2), -midpntvec(:,1)]));
                  
                  % for the bottom segment nodes, at the boundary of the region, we
                  % need to find the intercept of the sides of the geometry and a radial
                  % line going from the origin
                  xy1 = [ transnodes(segaboveind,:), transnodes(segaboveind,:) - 10*midpntvec; ...
                          transnodes(segaboveind+1,:), transnodes(segaboveind+1,:) - 10*midpntvec; ...
                          0, 0, 0, 0;
                          leftnodes(2,:), rightnodes(2,:) ];
              
              end
                      
          elseif i == 2
          
              if rectcentres(end,2) < spacecentres(end,2)
          
                  spacesegbelowind = find (nodes(:,2) < spacecentres(end,2));
                  spacesegbelowind = spacesegbelowind(end-1);
              
                  % get mid points of the horizontal segments below 
                  midpoints = [ xmid, nodes(spacesegbelowind,2); ];
                  
                  % transform to lie on arc
                  [midpoints(:,1), midpoints(:,2)] = pol2cart(midpoints(:,2), midpoints(:,1));
                  
                  % get the midpoint of the other edge of the magnet by rotating
                  % the top mipoint about (0,0) by -thetamag radians (clockwise with
                  % this formulation)
                  rotM = [ cos(-thetaspace)  -sin(-thetaspace); 
                           sin(-thetaspace)   cos(-thetaspace) ];
                           
                  midpoints = [ midpoints * rotM;
                                midpoints; ] ;
              
                  % get vector between midpoints pointing from lower to upper
                  midpntvec = midpoints(1,:) - midpoints(2,:);
                  
                  % shift the upper midpoint down to lie at the split point of the
                  % magnet
                  midpoints(1,:) = midpoints (2,:) + midpntvec * ((2 * thetapole - nodes(spacesegbelowind,2)) / thetaspace);
                  
                  % for the top segment nodes, at the boundary of the region, we
                  % need to find the intercept of the sides of the geometry and a radial
                  % line going from the origin through the top of the geometry
                  xy1 = [ transnodes(spacesegbelowind,:), transnodes(spacesegbelowind,:) + 10*midpntvec; ...
                          transnodes(spacesegbelowind+1,:), transnodes(spacesegbelowind+1,:) + 10*midpntvec; ...
                          transnodes(spacesegbelowind,:), transnodes(spacesegbelowind,:) + 10*midpntvec];
                      
              else
                  % get mid points of the horizontal segments below 
                  midpoints = [ xmid, nodes(segbelowind,2); ];
                  
                  % transform to lie on arc
                  [midpoints(:,1), midpoints(:,2)] = pol2cart(midpoints(:,2), midpoints(:,1));
                  
                  % get the midpoint of the other edge of the magnet by rotating
                  % the top mipoint about (0,0) by -thetamag radians (clockwise with
                  % this formulation)
                  rotM = [ cos(-thetamag)  -sin(-thetamag); 
                           sin(-thetamag)   cos(-thetamag) ];
                           
                  midpoints = [ midpoints * rotM;
                                midpoints; ] ;
              
                  % get vector between midpoints pointing from lower to upper
                  midpntvec = midpoints(1,:) - midpoints(2,:);
                  
                  % shift the upper midpoint down to lie at the split point of the
                  % magnet
                  midpoints(1,:) = midpoints (2,:) + midpntvec * ((2 * thetapole - nodes(segbelowind,2)) / thetamag);
                  
                  % get new node locations of magnets on left and right
                  leftnodes = bsxfun(@plus, midpoints, rmag/2 * unit ([-midpntvec(:,2), midpntvec(:,1)]));
                  rightnodes = bsxfun(@plus, midpoints, rmag/2 * unit ([midpntvec(:,2), -midpntvec(:,1)]));
                  
                  % for the top segment nodes, at the boundary of the region, we
                  % need to find the intercept of the sides of the geometry and a radial
                  % line going from the origin through the top of the geometry
                  xy1 = [ transnodes(segbelowind,:), transnodes(segbelowind,:) + 10*midpntvec; ...
                          transnodes(segbelowind+1,:), transnodes(segbelowind+1,:) + 10*midpntvec; ...
                          transnodes(segaboveind,:), transnodes(segaboveind,:) + 10*midpntvec];              
              end
          
          end

      
      elseif size (rectcentres,1) == 3
          % one of the magnets is split. This makes finding the appropriate
          % vector direction more difficult as the existing geometry edge
          % segments will not be in the right direction
          
          if i == 1
              
              % get mid points of the horizontal segments above 
              midpoints = [ xmid, nodes(segaboveind,2); ];
              
              % transform to lie on arc
              [midpoints(:,1), midpoints(:,2)] = pol2cart(midpoints(:,2), midpoints(:,1));
              
              % get the midpoint of the other edge of the magnet by rotating
              % the top mipoint about (0,0) by thetamag radians (anticlockwise
              % with this formulation)
              rotM = [ cos(thetamag)  -sin(thetamag); 
                       sin(thetamag)   cos(thetamag) ];
                       
              midpoints = [ midpoints; 
                            midpoints * rotM ];
                            
              % get vector between midpoints pointing from lower to upper
              midpntvec = midpoints(1,:) - midpoints(2,:);
              
              % shift the lower midpoint up to lie at the split point of the magnet
              midpoints(2,:) = midpoints (1,:) - midpntvec * (nodes(segaboveind,2) / thetamag);
              
              % get new node locations of magnets on left and right
              leftnodes = bsxfun(@plus, midpoints, rmag/2 * unit ([-midpntvec(:,2), midpntvec(:,1)]));
              rightnodes = bsxfun(@plus, midpoints, rmag/2 * unit ([midpntvec(:,2), -midpntvec(:,1)]));
              
              % for the bottom segment nodes, at the boundary of the region, we
              % need to find the intercept of the sides of the geometry and a radial
              % line going from the origin
              xy1 = [ transnodes(segaboveind,:), transnodes(segaboveind,:) - 10*midpntvec; ...
                      transnodes(segaboveind+1,:), transnodes(segaboveind+1,:) - 10*midpntvec; ...
                      leftnodes(2,:), rightnodes(2,:) ];

          elseif i == 2
              % can use normal method
              
              % get mid points of horizontal segments above and below
              midpoints = [ xmid, nodes(segaboveind,2); ...
                            xmid, nodes(segbelowind,2); ];
              
              % transform to lie on arc
              [midpoints(:,1), midpoints(:,2)] = pol2cart(midpoints(:,2), midpoints(:,1));
              
              % get vector between midpoints pointing from lower to upper
              midpntvec = midpoints(1,:) - midpoints(2,:);
              
              % get new node locations of magnets on left and right
              leftnodes = bsxfun(@plus, midpoints, rmag/2 * unit ([-midpntvec(:,2), midpntvec(:,1)]));
              rightnodes = bsxfun(@plus, midpoints, rmag/2 * unit ([midpntvec(:,2), -midpntvec(:,1)]));
              
          elseif i == 3
              
              % get mid points of the horizontal segments below 
              midpoints = [ xmid, nodes(segbelowind,2); ];
              
              % transform to lie on arc
              [midpoints(:,1), midpoints(:,2)] = pol2cart(midpoints(:,2), midpoints(:,1));
              
              % get the midpoint of the other edge of the magnet by rotating
              % the top mipoint about (0,0) by -thetamag radians (clockwise with
              % this formulation)
              rotM = [ cos(-thetamag)  -sin(-thetamag); 
                       sin(-thetamag)   cos(-thetamag) ];
                       
              midpoints = [ midpoints * rotM;
                            midpoints; ] ;
          
              % get vector between midpoints pointing from lower to upper
              midpntvec = midpoints(1,:) - midpoints(2,:);
              
              % shift the upper midpoint down to lie at the split point of the
              % magnet
              midpoints(1,:) = midpoints (2,:) + midpntvec * ((2 * thetapole - nodes(segbelowind,2)) / thetamag);
              
              % get new node locations of magnets on left and right
              leftnodes = bsxfun(@plus, midpoints, rmag/2 * unit ([-midpntvec(:,2), midpntvec(:,1)]));
              rightnodes = bsxfun(@plus, midpoints, rmag/2 * unit ([midpntvec(:,2), -midpntvec(:,1)]));
              
              % for the top segment nodes, at the boundary of the region, we
              % need to find the intercept of the sides of the geometry and a radial
              % line going from the origin through the top of the geometry
              xy1 = [ transnodes(segbelowind,:), transnodes(segbelowind,:) + 10*midpntvec; ...
                      transnodes(segbelowind+1,:), transnodes(segbelowind+1,:) + 10*midpntvec; ...
                      transnodes(segaboveind,:), transnodes(segaboveind,:) + 10*midpntvec; ...
                      leftnodes(2,:), rightnodes(2,:)];
              
          end
      
      end
      
      % insert into the translated nodes matrix
      transnodes([segaboveind,segbelowind],:) = leftnodes;
      transnodes([segaboveind+1,segbelowind+1],:) = rightnodes;
      
      if i == 1
      
          % for the bottom segment nodes, at the boundary of the region, we
          % need to find the intercept of the sides of the geometry and a radial
          % line going from the origin
          
%          xy1 = [ transnodes(segaboveind,:), transnodes(segaboveind,:) - 10*midpntvec; ...
%                  transnodes(segaboveind+1,:), transnodes(segaboveind+1,:) - 10*midpntvec; ];
              
          xy2 = repmat([ 0, 0, 1000, 0;], size(xy1,1), 1);
                  
          intersections = lineSegmentIntersect (xy1,xy2);
          
          if transnodes(segbelowind,2) < Inputs.Tol
              % a magnet node is poking out the bottom, to a distance greater 
              % than the allowed tolerance, this requires modification of the
              % geometry at both top and bottom. 
              
              % an additional two nodes and two links must be added at the bottom
              % right, and an additional node and two links must be added at the
              % top left
              
              origtnode = [ transnodes(1,:);
                            intersections.intMatrixX(4,4), intersections.intMatrixY(4,4); ];
                            
              rotM = [ cos(-2*thetapole)  -sin(-2*thetapole); 
                       sin(-2*thetapole)   cos(-2*thetapole) ];
                       
              newtnode = origtnode * rotM;
              
              transnodes = [ intersections.intMatrixX(1,1), intersections.intMatrixY(1,1); ...
                             intersections.intMatrixX(2,2), intersections.intMatrixY(2,2); ...
                             intersections.intMatrixX(4,4), intersections.intMatrixY(4,4); ...
                             transnodes(2:end-2,:); ...
                             newtnode(1,:); ...
                             (roffset + rmag/2) * [cos(2*thetapole), sin(2*thetapole)]; ...
                             (roffset - rmag/2) * [cos(2*thetapole), sin(2*thetapole)]; ...
                             newtnode(2,:) ];

              ntransnodes = size(transnodes,1);
              
              % increment all the link node ids
              links(links > elcount.NNodes) = links(links > elcount.NNodes) + 2;
              
              % add new links
              links = [ links(1,1), elcount.NNodes + 2; ...
                        elcount.NNodes + 2, elcount.NNodes + 1; ...
                        elcount.NNodes + 1, elcount.NNodes + 3; ...
                        elcount.NNodes + 2, elcount.NNodes + 3; ...
                        links(2:end-3,:); ...
                        min(links(end-2,:)), elcount.NNodes+ntransnodes-4; ...
                        min(links(end-1,:)), elcount.NNodes+ntransnodes-3; ...
                        elcount.NNodes+ntransnodes-4, elcount.NNodes+ntransnodes-1;
                        elcount.NNodes+ntransnodes-4, elcount.NNodes+ntransnodes-2;
                        elcount.NNodes+ntransnodes-1, elcount.NNodes+ntransnodes-3; 
                        elcount.NNodes+ntransnodes-2, elcount.NNodes+ntransnodes-1; ];
                        
              elcount.NNodes = elcount.NNodes + 2;
                             
              % insert placeholder node of nans in original nodes matrix
              nodes = [ nodes(1,:);
                        nan, nan;
                        nan, nan;
                        nodes(2:end,:) ];
                        
              % create new label location
              extraspacecentres = tricenter ( [transnodes(2,:), 0], ...
                                              [transnodes(3,:), 0], ...
                                              [transnodes(4,:), 0] ).';
                                 
              % set the flag indicating the modification to true
              edgemods = true;
              
          elseif transnodes(segbelowind,2) < 0
              % a magnet node is poking out the bottom, to a distance less 
              % than the allowed tolerance, so just shift it to the intercept
              transnodes(1,:) = [intersections.intMatrixX(1,1), intersections.intMatrixY(1,1)];
              transnodes(2,:) = [intersections.intMatrixX(2,2), intersections.intMatrixY(2,2)];
          else
              % for the bottom segment nodes, at the boundary of the region, we
              % need to find the intercept of the sides of the geomtry and a
              % radial line going from the origin, and move the nodes to these
              % points
              
              % set the new node locations for the bottom segment
              transnodes(1,:) = [ intersections.intMatrixX(1,1), intersections.intMatrixY(1,1) ];
              transnodes(2,:) = [ intersections.intMatrixX(2,2), intersections.intMatrixY(2,2) ];
              
          end
          
      elseif i == size (rectcentres,1)
                  
          xy2 = repmat([ 0, 0, 10 * [cos(2*thetapole), sin(2*thetapole)]], 3, 1);
                      
          intersections = lineSegmentIntersect (xy1,xy2);
          
          if intersections.intAdjacencyMatrix(3,3) == 1 && intersections.intNormalizedDistance1To2(3,3) > Inputs.Tol
              % a magnet node is poking out the top, to a distance greater 
              % than the allowed tolerance, this requires modification of the
              % geometry at both top and bottom. 
              
              % an additional node and two links must be added at the bottom
              % right, and an additional node and two links must be added at the
              % top left
              
              transnodes = [ transnodes(1:end-2,:); ...
                             intersections.intMatrixX(2,2), intersections.intMatrixY(2,2); ...
                             intersections.intMatrixX(4,4), intersections.intMatrixY(4,4); ...
                             intersections.intMatrixX(1,1), intersections.intMatrixY(1,1); ...
                             transnodes(end,:); ];
                             
              % increment all the link node ids
%              links(links > nnodes) = links(links > nnodes) + 2;

              % add new link
              links = [ links(1:end-1,:);
                        nnodes + 2, nnodes + 1; ...
                        nnodes + 2, nnodes + 1; ...
                        links(end,1)+2, nnodes + 2 ];
                        
              % insert placeholder node of nans
              nodes = [ nodes(1,:);
                        nan, nan;
                        nan, nan;
                        nodes(2:end,:) ];
              
          elseif intersections,intAdjacencyMatrix(3,3) == 1 && intersections.intNormalizedDistance1To2(3,3) > 0
              % a magnet node is poking out the top, to a distance less 
              % than the allowed tolerance, so just shift it to the intercept
              transnodes(end-1,:) = [intersections.intMatrixX(1,1), intersections.intMatrixY(1,1)];
              transnodes(end,:) = [intersections.intMatrixX(2,2), intersections.intMatrixY(2,2)];
          else
              % for the top segment nodes, at the boundary of the region, we
              % need to find the intercept of the sides of the geomtry and a
              % radial line going from the origin, and move the nodes to these
              % points
              
              % set the new node locations for the bottom segment
              transnodes(end-1,:) = [ intersections.intMatrixX(1,1), intersections.intMatrixY(1,1) ];
              transnodes(end,:) = [ intersections.intMatrixX(2,2), intersections.intMatrixY(2,2) ];
              
          end
      
      end
      
    end
    
    
%   
%    nodes(1:2:end-1,:) = leftnodes;
%    nodes(2:2:end,:) = leftnodes; 
%    
%    % for the top and bottom segment nodes, at the boundary of the region, we
%    % need to find the intercept of the sides of the magnet and a radial line
%    % going from the origin, and move the nodes to these points
%    xy1 = [ nodes(3,:), nodes(3,:) - 3*midpntvecs(1,:); ...
%            nodes(4,:), nodes(4,:) - 3*midpntvecs(1,:); ...
%            nodes(end-1,:), nodes(end-1,:) + 3*midpntvecs(end,:); ...
%            nodes(end,:), nodes(end,:) + 3*midpntvecs(end,:); ];
%            
%    xy2 = [ 0, 0, 10 * midpoints(1,1), 0; ...
%            0, 0, 10 * midpoints(1,1), 0; ...
%            0, 0, 10 * midpoints(end,:); ...
%            0, 0, 10 * midpoints(end,:); ];
%            
%    intersections = lineSegmentIntersect (xy1,xy2);
%    
%    % get the new node locations
%    nodes(1,:) = [ intersections.intMatrixX(1,1), intersections.intMatrixY(1,1) ];
%    nodes(2,:) = [ intersections.intMatrixX(2,2), intersections.intMatrixY(2,2) ];
%    nodes(end-1,:) = [ intersections.intMatrixX(3,3), intersections.intMatrixY(3,3) ];
%    nodes(end,:) = [ intersections.intMatrixX(4,4), intersections.intMatrixY(4,4) ];
    
    % add all the nodes to the problem
    for i = 1:size(transnodes,1)
        FemmProblem = addnodes_mfemm(FemmProblem, transnodes(i,1), transnodes(i,2), 'InGroup', Inputs.MagnetGroup);
    end
    
    % Periodic boundary at bottom
    FemmProblem = addsegments_mfemm(FemmProblem, links(1,1), links(1,2), ...
        'BoundaryMarker', FemmProblem.BoundaryProps(boundind).Name, ...
        'InGroup', Inputs.MagnetGroup);
        
    % Add all the segments except the top and bottom which will have
    % boundary properties we must set manually
    FemmProblem = addsegments_mfemm(FemmProblem, links(2:end-1,1), links(2:end-1,2), ...
        'InGroup', Inputs.MagnetGroup);

    % Periodic boundary at top
    FemmProblem = addsegments_mfemm(FemmProblem, links(end,1), links(end,2), ...
        'BoundaryMarker', FemmProblem.BoundaryProps(boundind).Name, ...
        'InGroup', Inputs.MagnetGroup);
    
    % Now add the labels
    [rectcentres(:,1), rectcentres(:,2)] = pol2cart(rectcentres(:,2), rectcentres(:,1));
    [spacecentres(:,1), spacecentres(:,2)] = pol2cart(spacecentres(:,2), spacecentres(:,1));
    r = magn (rectcentres(1,1:2));
    
    if ~isempty(Inputs.MagnetMaterial)
        
        magblockinds = [];
        
        % add the magnet labels
        for i = 1:size(rectcentres, 1)

            labelpnt = unit (rectcentres(i,1:2)) * sqrt (r^2 .* (cos (thetamag/2))^2);

            if rectcentres(i,3)

                [FemmProblem, magblockinds(end+1,1)] = addblocklabel_mfemm(FemmProblem, labelpnt(1), labelpnt(2), ...
                                                'BlockType', FemmProblem.Materials(Inputs.MagnetMaterial).Name, ...
                                                'MaxArea', Inputs.MeshSize, ...
                                                'MagDir', Inputs.MagDirections{2}, ...
                                                'InGroup', Inputs.MagnetGroup);
                                                
                if edgemods
                
                end
                
            else

                [FemmProblem, magblockinds(end+1,1)] = addblocklabel_mfemm(FemmProblem, labelpnt(1), labelpnt(2), ...
                                                'BlockType', FemmProblem.Materials(Inputs.MagnetMaterial).Name, ...
                                                'MaxArea', Inputs.MeshSize, ...
                                                'MagDir', Inputs.MagDirections{1}, ...
                                                'InGroup', Inputs.MagnetGroup);

                if edgemods
                
                end
                
            end
            
            magblockinds(end, 2) = rectcentres(i,3);

        end
    
    end
    
    if ~isempty(Inputs.SpaceMaterial)
        
        % Add the other lables
        for i = 1:size(spacecentres, 1)
            
            labelpnt = unit (spacecentres(i,:)) * sqrt (r^2 .* (cos ((thetapole - thetamag)/2))^2);

            FemmProblem = addblocklabel_mfemm(FemmProblem, labelpnt(1), labelpnt(2), ...
                                            'BlockType', FemmProblem.Materials(Inputs.SpaceMaterial).Name, ...
                                            'MaxArea', Inputs.MeshSize, ...
                                            'InGroup', Inputs.SpaceGroup);

        end
        
        if edgemods
        
            FemmProblem = addblocklabel_mfemm(FemmProblem, extraspacecentres(1), extraspacecentres(2), ...
                                            'BlockType', FemmProblem.Materials(Inputs.SpaceMaterial).Name, ...
                                            'MaxArea', Inputs.MeshSize, ...
                                            'InGroup', Inputs.SpaceGroup);
        
        end
    
    end

    
end