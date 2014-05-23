function [nodes, links, cornernodes, shoegaplabelloc, coillabelloc, vertlinkinds, toothlinkinds] ...
    = internalslotnodelinks(ycoil, yshoegap, xcore, xcoil, xshoebase, xshoegap, ylayers, tol, varargin)
% creates a set of node, links and label locations for the points making up
% the inside of a coil slot in an iron cored machine
%
% Syntax
%
% [nodes, links, cornernodes, shoegaplabelloc, coillabelloc, vertlinkinds, toothlinkinds] = ...
%        internalslotnodelinks(ycoil, yshoegap, xcore, xcoil, xshoebase, xshoegap, layers, tol)
%
% Description
%
% internalslotnodelinks is a low-level function used to create a set of
% node, links and label locations for the points making up the internal
% surface of a single coil slot in an iron cored machine. It is intended
% that these slot locations are replicated and joined by a higher level
% function to create a machine stator component. The created slot has a
% shape similar to that shown in Figure 1.
%
%                   xcore    xcoil   xshoebase
%                   <----><----------><-->
%   ^               |     |              _ 
%                   |     |             / |
%                   |     |___________ /  |
%                   |                     |
%                   |                     |
%                   |                     |
%                   |      ___________    |
%                ^  |     |            \  |
%                :  |     |             \_|
%                :  |     |                ^        
%         ycoil  :  |     |                : yshoegap                
%   ^            :  |     |              _ v
%                :  |     |             / |
%                v  |     |___________ /  |
%                   |                     |
%                   |                     |
%                   |                     |
%                   |      ___________    |
%                   |     |            \  |
%                   |     |             \_|
%                                       <->
%                                      xshoegap
%
%                    Figure 1 - slot shapes
%         
%
% The actual output of internalslotnodelinks is shown in Figure 2, which is
% the internal surface of the slot.
%
%                                                \ / outer corner node
%                  |                             / \
%                   
%                  |
%                   xcore    xcoil   xshoebase
%                  |<----><----------><-->
%                                        \ / cn3
%                  |       cn4           /|\
%                        \ /__________    |
%                ^ |     /|\           \  |
%                :        |             \_|
%                : |      |                ^        
% - - - - ycoil  : - - - -|- - - - - - - - : yshoegap - - - - x = 0  
%                : |      |              _ v
%                :        |             / |
%                v |     \|/__________ /  |
%                        / \              |
%                  |      cn1            \|/
%                                        / \ cn2
%                  |                    <->
%                                      xshoegap
%                y = 0
%
%       Figure 2 - Acutal output of internalslotnodelinks
%    
% internalslotnodelinks produces a generic slot shape in 2D cartesian
% coordinates with the centreline of the slot along the line x = 0 and the
% back of the slot yoke along the line y = 0 as in Fig. 2. This slot may
% require modification (e.g. rotation, or mirroring) by higher level
% routines for particular machine application, e.g. a radial flux machine.
%
% The slot can also be split into multiple layers in the x direction. Each
% layer will be separated by a link, and a label location will be provided
% for each layer.
%
% Input
%
%  ycoil - internal slot width in y direction
% 
%  yshoegap - distance between shoe tips in y direction. This can be any
%    dimension from zero up to 
% 
%  xcore - thickness of slot yoke in x direction (in practice this is the
%    distance of the base of the drawn slot shape from the y-axis) 
% 
%  xcoil - internale slot height in y direction, if layers are specified
%    it is this dimension which is split into layers.
% 
%  xshoebase - thichness of the shoe at the point it meets the tooth yoke 
% 
%  xshoegap - thickness of the shoe at its tip, this can be zero, resulting
%    in a pointed tip. Note that this can cause meshing issues depending on
%    the angle of this pointed tip.
% 
%  layers - the number of layers of coils in the slots. The slot will be
%    split up into  the specified number of layers in the x direction.
% 
%  tol - tolerance determining whether certain dimensions are joined or
%    not, e.g. whether yshoegap should be treated as zero. Default is 1e-5.
%
%
% Output
%
% nodes - (n x 2) matrix of x and y coordinates of the nodes making up the
%  internal dimensions of the slots
%
% links - (n x 2) matrix of zero-based links between the nodes. The values
%   in links are the rows of 'nodes' counting from zero. 
% 
% cornernodes - vector containin the zero-based rows of the nodes at the
%   outer corners of the slot shape. This are supplied clockwise from the
%   bottom left, as shown in Figure 2, such that the vector contains
%
%   [ cn1, cn2, cn3, cn4 ]
%  
% shoegaplabelloc - (n x 2) matrix of label locations of any gaps created
%   by the shoe geometry. If yshoegap (and yshoebase) is greater than zero,
%   there will be an air-filled rectangular space between the shoe tips. A
%   suitable label location for this region is provided here. If the shoes
%   have a wider base than tip, there will be an additional trianglular or
%   trapezoidal gap between the top of the coil and the start of the space
%   between the shoe tips, and another suitable label location will be
%   provided. If xshoebase is zero (so there are no shoes), shoegaplabelloc
%   will be empty.
% 
% coillabelloc - (n x 2) matrix of coil label locations. The number of rows
%   will depend on the value of layers, with one label for every layer.
%
% vertlinkinds - vector of (one-based) indices of the links which should be 
%   converted to arcs when using internalslotnodelinks to make a radial 
%   design, e.g the outer tooth surface links.
%
% toothlinkinds - vector of indices of the links that make up the internal
%   surface of the tooth, not including the links creating divisions in the
%   cross-section of the slot of coils, or the marking the air region
%   between the soes if there is a gap.
%

% Created by Richard Crozier 2012

    if nargin == 0 && nargout == 1
        % return function handles to subfunctions for testing purposes
        nodes = {@shoecurvepoints, @basecurvepoints};
        return;
    end
    
    options.CoilBaseFraction = 0.05;
    options.ShoeCurveControlFrac = 0.5;
    options.SplitX = false;
    
    options = parseoptions (options, varargin);
    
    if xcoil < 5*tol
        error ('xcoil cannot be less than five times the tolerance.')
    end
    
    if options.CoilBaseFraction > 1 || options.CoilBaseFraction < 0
        error ('optional CoilBaseFraction must lie between 0 and 1.')
    end
    
    xcoilbase = options.CoilBaseFraction * xcoil;
    
    if xcoilbase < 2*tol
        xcoilbase = 2*tol;
    end
    
    xcoilbody = xcoil - xcoilbase;
    
    % there are no shoes so return an empty label
    shoegaplabelloc = [];
    
    if yshoegap > ycoil
        yshoegap = ycoil;
    end
    
    if numel (ycoil) == 1
        ycoilshoe = ycoil;
        ycoilbase = ycoil;
    elseif numel (ycoil) == 2
        ycoilshoe = ycoil(1);
        ycoilbase = ycoil(2);
    else
        error ('ycoil must contain only one or two elements');
    end
        
    if xshoebase < tol
        % there is no shoe in this case, so just make a line for the
        % coil flush with the tooth surface
        
        nodes = [ xcore + xcoilbase + xcoilbody, ycoilshoe/2;
                  xcore + xcoilbase + xcoilbody, -ycoilshoe/2;
                  xcore + xcoilbase + xcoilbody, 0 ];
        
        links = [ 1, 2;
                  2, 0 ];
              
        shoex = [];
        shoey = [];
        topshoenids = [];
        botshoenids = [];

        topinnershoenode = 0;
        botinnershoenode = 1;
        
        topouternode = 0;
        botouternode = 1;
        
        midshoenode = 2;
        
        vertlinkinds = [1, 2];
        
        toothlinkinds = [];
                
    else
        % there is a shoe on top of the tooth

        if yshoegap < tol
            % The shoe is joined in the middle seems odd, but we'll allow
            % for it here, who knows what applications the future might
            % have for such a geometry?

            if xshoegap < tol
                
                % the tip of the shoe meets in a sharp point
                %
                % In this case three main nodes are required to draw the
                % shoe
 
                % calculate intermediate points between the shoe and top of
                % coil region which create a curved surface
                [shoex, shoey, shoeQx, shoeQy] = shoecurvepoints ( options.ShoeCurveControlFrac, ...
                                                   xcore, ...
                                                   xcoil, ...
                                                   xshoebase, ...
                                                   xshoegap, ...
                                                   ycoilbase, ...
                                                   ycoilshoe, ...
                                                   yshoegap );
                                 
                ncurvepnts = numel(shoex);
                
                % the top shoe nodes
                nodes = [ xcore + xcoilbase + xcoilbody, ycoilshoe/2; ...
                          xcore + xcoilbase + xcoilbody + xshoebase, ycoilshoe/2; ...
                          xcore + xcoilbase + xcoilbody + xshoebase, 0; ...
                          xcore + xcoilbase + xcoilbody, -ycoilshoe/2; ...
                          xcore + xcoilbase + xcoilbody + xshoebase, -ycoilshoe/2; ...
                          shoex', shoey'; ...
                          shoex', -shoey' ];                       
                
                topshoenids = 5:4+ncurvepnts;
                botshoenids = 5+ncurvepnts:4+2*ncurvepnts;
                
                links = [ 2, 1; ...
                          0, topshoenids(1); ...
                          [topshoenids(1:end-1)', topshoenids(2:end)' ]; ... % join the curve nodes
                          topshoenids(end), 2; ...
                          3, botshoenids(1); ...
                          [botshoenids(1:end-1)', botshoenids(2:end)' ]; ...
                          botshoenids(end), 2; ...
                          2, 4 ];

                topinnershoenode = 0;
                botinnershoenode = 3; 
                
                topouternode = 1;
                botouternode = 4;
                
                midshoenode = 2;
                
                vertlinkinds = [ 1, size(links,1)];

            else
                % there is a blunt edge on the shoe
                %
                % In this case 4 nodes are required to draw the shoe
                % corners
                
                % calculate intermediate points between the shoe and top of
                % coil region which create a curved surface
                [shoex, shoey, shoeQx, shoeQy] = shoecurvepoints ( options.ShoeCurveControlFrac, ...
                                                   xcore, ...
                                                   xcoil, ...
                                                   xshoebase, ...
                                                   xshoegap, ...
                                                   ycoilbase, ...
                                                   ycoilshoe, ...
                                                   yshoegap );

                ncurvepnts = numel(shoex);
                
                nodes = [ xcore + xcoilbase + xcoilbody, ycoilshoe/2; ...
                          xcore + xcoilbase + xcoilbody + xshoebase, ycoilshoe/2; ...
                          xcore + xcoilbase + xcoilbody + xshoebase, 0; ...
                          xcore + xcoilbase + xcoilbody + xshoebase - xshoegap, 0];

                % now draw the other shoe 
                nodes = [ nodes; 
                          nodes(2,1), -nodes(2,2); ...
                          nodes(1,1), -nodes(1,2); ...
                          shoex', shoey';
                          shoex', -shoey'];

                topshoenids = 6:5+ncurvepnts;
                botshoenids = 6+ncurvepnts:5+2*ncurvepnts;
                
                links = [ 2, 1;
                          2, 3;
                          0, topshoenids(1); 
                          [topshoenids(1:end-1)', topshoenids(2:end)' ]; ... % join the curve nodes
                          topshoenids(end), 3; ...
                          5, botshoenids(1); ...
                          [botshoenids(1:end-1)', botshoenids(2:end)' ]; ...
                          botshoenids(end), 3;
                          4, 2];

                topinnershoenode = 0;
                botinnershoenode = 5; 
                
                topouternode = 1;
                botouternode = 4;
                
                midshoenode = 3;
                
                vertlinkinds = [ 1, size(links,1)];
                
                
            end
            
            % all the links are the tooth surface at this stage
            toothlinkinds = links;

        else
            % we must acount for the space between shoes in the model 
            if xshoegap < tol
                % the tip of the shoe ends in a sharp point

                [shoex, shoey, shoeQx, shoeQy] = shoecurvepoints ( options.ShoeCurveControlFrac, ...
                                                                   xcore, ...
                                                                   xcoil, ...
                                                                   xshoebase, ...
                                                                   xshoegap, ...
                                                                   ycoilbase, ...
                                                                   ycoilshoe, ...
                                                                   yshoegap );
                                            
                ncurvepnts = numel(shoex);
                
                % the top shoe nodes
                nodes = [ xcore + xcoilbase + xcoilbody, ycoilshoe/2; ...
                          xcore + xcoilbase + xcoilbody + xshoebase, ycoilshoe/2; ...
                          xcore + xcoilbase + xcoilbody + xshoebase, yshoegap/2; ...
                          xcore + xcoilbase + xcoilbody, -ycoilshoe/2; ...
                          xcore + xcoilbase + xcoilbody + xshoebase, -ycoilshoe/2; ...
                          xcore + xcoilbase + xcoilbody + xshoebase, -yshoegap/2; ...
                          xcore + xcoilbase + xcoilbody + xshoebase, 0;
                          shoex', shoey'; ...
                          shoex', -shoey' ];
                      
                topshoenids = 7:6+ncurvepnts;
                botshoenids = 7+ncurvepnts:6+2*ncurvepnts;
                
                links = [ 1, 2;
                          0, topshoenids(1); 
                          [topshoenids(1:end-1)', topshoenids(2:end)' ]; ... % join the curve nodes
                          topshoenids(end), 2; 
                          3, botshoenids(1); 
                          [botshoenids(1:end-1)', botshoenids(2:end)' ]; 
                          botshoenids(end), 5; 
                          4, 5;
                          5, 6
                          6, 2 ];

                topinnershoenode = 0;
                botinnershoenode = 3; 
                
                topouternode = 1;
                botouternode = 4;
                
                midshoenode = 6;
                
                vertlinkinds = [ 1, size(links,1)-2];
                
                % tooth links are all links except the links closing the
                % air gap
                toothlinkinds = 1:size(links,1);
                toothlinkinds([size(links,1), size(links,1)-1]) = [];
                          
            else
                % blunt edged shoes with gap
                
                % calculate intermediate points between the shoe and top of
                % coil region which create a curved surface
                [shoex, shoey, shoeQx, shoeQy] = shoecurvepoints ( options.ShoeCurveControlFrac, ...
                                                   xcore, ...
                                                   xcoil, ...
                                                   xshoebase, ...
                                                   xshoegap, ...
                                                   ycoilbase, ...
                                                   ycoilshoe, ...
                                                   yshoegap );

                ncurvepnts = numel(shoex);
                
                % there is a blunt edge on the shoe
                nodes = [ xcore + xcoilbase + xcoilbody, ycoilshoe/2; ...
                          xcore + xcoilbase + xcoilbody + xshoebase, ycoilshoe/2; ...
                          xcore + xcoilbase + xcoilbody + xshoebase, yshoegap/2; ...
                          xcore + xcoilbase + xcoilbody + xshoebase - xshoegap, yshoegap/2];
                      
                % now draw the other shoe 
                nodes = [ nodes; 
                          nodes(:,1), -nodes(:,2); 
                          xcore + xcoilbase + xcoilbody + xshoebase - xshoegap, 0;
                          shoex', shoey';
                          shoex', -shoey' ];    
                      
                topshoenids = 9:8+ncurvepnts;
                botshoenids = 9+ncurvepnts:8+2*ncurvepnts;
                
                links = [ 2, 1;
                          5, 6;
                          2, 3;
                          0, topshoenids(1);
                          [topshoenids(1:end-1)', topshoenids(2:end)'  ]; ... % join the curve nodes
                          topshoenids(end), 3;
                          6, 7;
                          4, botshoenids(1);
                          [botshoenids(1:end-1)', botshoenids(2:end)'  ]; ... % join the curve nodes
                          botshoenids(end), 7 ];

                % tooth links are all links except the links closing the
                % air gap
                toothlinkinds = 1:size(links,1);
                
                % make the shoe gap region
                links = [ links; 
                          6, 2; 
                          7, 8;
                          8, 3 ];
    
                % get the centre of the shoe gap
                shoegaplabelloc = rectcentre(nodes(8,:), nodes(3,:));
                       
                topinnershoenode = 0;
                botinnershoenode = 4; 
                
                topouternode = 1;
                botouternode = 5;
                
                midshoenode = 8;
                
                vertlinkinds = [ 1, 2, size(links,1)-2 ];
                
            end

        end

    end
    
    % make the curve at the base of the slot
    [basex, basey, baseQx, baseQy, m, c] = basecurvepoints (xcore, xcoilbase, xcoilbody, ycoilbase, ycoilshoe);
    
    ncurvepnts = numel (basex);
    
    % this will be the id of the next node we add, the (xcore, 0) node
    startbasenodeid = size (nodes, 1);
    
    nodes = [ nodes;
              xcore, 0; 
              basex', basey'; 
              basex', -basey';
              xcore+xcoilbase, ycoilbase/2;
              xcore+xcoilbase, -ycoilbase/2 ];
          
    endbasenodeid = size (nodes, 1) - 1;

    linksize = size(links,1);
    
    links = [ links;
              startbasenodeid, startbasenodeid+1; 
              [(startbasenodeid+1:startbasenodeid+ncurvepnts-1)', (startbasenodeid+2:startbasenodeid+ncurvepnts)']
              startbasenodeid+ncurvepnts, endbasenodeid - 1;
              startbasenodeid, startbasenodeid+ncurvepnts+1;
              [(startbasenodeid+ncurvepnts+1:startbasenodeid+2*ncurvepnts-1)', (startbasenodeid+ncurvepnts+2:startbasenodeid+2*ncurvepnts)'] 
              startbasenodeid+2*ncurvepnts, endbasenodeid  ];
    
    % store the tooth surface link indices
    toothlinkinds = [toothlinkinds, linksize+1:size(links,1)];
    
    % get the starting inner nodes for making links
    lastinnernodes = [ topinnershoenode, botinnershoenode ];
    
    % starting position for layers, the top of the outer layer
    
    % calculate area of various coil regions, we will operate on one half
    % of the winding to simplify finding the area under curves etc
    basex = [xcore, basex, xcore+xcoilbase];
    basey = [0, basey, ycoilbase/2];
    topbasenids = [(startbasenodeid:startbasenodeid+ncurvepnts), endbasenodeid - 1];
	botbasenids = [startbasenodeid, (startbasenodeid+ncurvepnts+1:startbasenodeid+2*ncurvepnts), endbasenodeid];
    
    basearea = trapz (basex, basey);
    bodyarea = trapzarea ( ycoilbase, ycoilshoe, xcoilbody ) / 2;
    shoearea = trapz ([ xcore+xcoil, shoex, xcore+xcoil+xshoebase-xshoegap ], ...
                      [ ycoilbase/2, shoey, yshoegap ] );
    
    totalarea = basearea + bodyarea + shoearea;
    
    % calculate the area of one winding layer, given the number of layers
    windingarea = totalarea / ylayers;
    
    if ylayers == 1
        % simple case, there is only one layer per slot, winding fills the
        % entire slot, but with the option of splitting into two layers in
        % the alternate direction so coils lie side-by-side one another in
        % the slot instead of on top of one another
        
        if options.SplitX
            
            % add a link creating the sides and a horizontal split
            links = [ links; 
                      startbasenodeid, midshoenode;
                      lastinnernodes(1), endbasenodeid - 1;
                      lastinnernodes(2),  endbasenodeid ];
            
            % get a suitible position for the coil labels
            coillabelloc = [ xcore + xcoilbase + xcoilbody/2, ycoilshoe/4;
                             xcore + xcoilbase + xcoilbody/2, -ycoilshoe/4 ];
        else
            % just add links making the sides
            links = [ links; 
                      lastinnernodes(1), endbasenodeid - 1;
                      lastinnernodes(2),  endbasenodeid ];
                  
            % get a suitible position for the coil label         
            coillabelloc = [ xcore + (xcoilbase + xcoilbody)/2, 0 ];
        end
        
        % update the tooth links
        toothlinkinds = [toothlinkinds, size(links, 1) - 1, size(links, 1)];
        
        cornernodes = [ lastinnernodes(2), botouternode, topouternode, lastinnernodes(1) ];
        
    else
        
        layersmade = 0;
        
        coillabelloc = [];
        
        basemidnode = startbasenodeid;
        
        layer_area_available = 0;
        
        lastlayerstartx = xcore;
        
        if windingarea < basearea
            % gobble up the base area, creating layers as necessary
            starttrapzind = 1;
                    
            for ind = 2:numel (basex)
                gobblearea = trapz (basex(starttrapzind:ind), basey(starttrapzind:ind));
                if gobblearea >= windingarea
                    % create a layer link and coil label location
                    
                    links = [ links; topbasenids(ind), botbasenids(ind) ];
                    
                    coillabelloc = [ coillabelloc; lastlayerstartx + (basex(ind) - lastlayerstartx)/2, 0 ];
                    layersmade = layersmade + 1;
                    
                    layer_area_available = basearea - gobblearea;
                    starttrapzind = ind;
                    lastlayerstartx = basex(ind);
                else
                    layer_area_available = gobblearea;
                end
            end
        else
            layer_area_available = basearea;
        end
        
        if (layersmade >= ylayers)
            % return if we've made enough layers, partly to sidestep
            % numerical issues
            cornernodes = [ lastinnernodes(2), botouternode, topouternode, lastinnernodes(1) ];
            return;
        end

        if windingarea < (layer_area_available + bodyarea)
            % we calculate the required position for x by finding the roots
            % of a quadratic obtained from the integral of the staight line
            % making up the straight side of the coil
            x1 = xcore + xcoilbase;
            
            lastlayernodeids = [topbasenids(end), botbasenids(end)];
            
            while 1
                
                bodyarealeft = trapzarea ( mxplusc (m, c, x1)*2, ...
                                           mxplusc (m, c, xcore + xcoil)*2, ...
                                           (xcore + xcoil) - x1 ) / 2;

                bodyarealeft = bodyarealeft + layer_area_available;
                
                if (windingarea > bodyarealeft) || (layersmade >= ylayers) ...

                    % link up the sides
                    links = [ links;
                              lastlayernodeids(1), topinnershoenode;
                              lastlayernodeids(2), botinnershoenode; ];
                    
                    % update the tooth links
                    toothlinkinds = [toothlinkinds, size(links, 1) - 1, size(links, 1)];
                    
                    % store the remaining area in the body available to
                    % a layer
                    layer_area_available = bodyarealeft;
                    
                    % break out of the while loop
                    break;
                    
                elseif (xshoebase < tol) && (layersmade == (ylayers - 1))
                    
                    % link up the sides
                    links = [ links;
                              lastlayernodeids(1), topinnershoenode;
                              lastlayernodeids(2), botinnershoenode; ];

                    % update the tooth links
                    toothlinkinds = [toothlinkinds, size(links, 1) - 1, size(links, 1)];

                    % put the remaining layer in the space between the last area
                    % boundary and the top
                    coillabelloc = [ coillabelloc; lastlayerstartx + ((xcore+xcoilbody+xcoilbase) - lastlayerstartx)/2, 0 ];

                    % return as there is no shoe area to gobble
                    cornernodes = [ lastinnernodes(2), botouternode, topouternode, lastinnernodes(1) ];
                    return;
                    
                else
                    % calculate the area of the body required to make up the
                    % remaining required winding area
                    At = windingarea - layer_area_available;

                    % find the roots of the polynomial to get the position of x
                    arearoots = roots ([ m/2, c, -(At + (m/2)*x1 + c*x1)]);

                    x2 = arearoots(arearoots > 0);

                    if isempty (x2)
                        error ('Impossible coil shape');
                    end

                    % add new nodes and link
                    newy = mxplusc (m, c, x2);

                    nodes = [ nodes; 
                              x2, newy; 
                              x2, -newy ];

                    thislayernodeids = [ size(nodes,1)-2, size(nodes,1)-1 ];

                    links = [ links;
                              lastlayernodeids(1), thislayernodeids(1);
                              lastlayernodeids(2), thislayernodeids(2);
                              thislayernodeids(1), thislayernodeids(2) ];
                          
                    % update the tooth links
                    toothlinkinds = [toothlinkinds, size(links, 1) - 2, size(links, 1) - 1];

                    coillabelloc = [ coillabelloc; lastlayerstartx + (x2-lastlayerstartx)/2, 0 ];
                    layersmade = layersmade + 1;
                    
                    lastlayernodeids = thislayernodeids;

                    lastlayerstartx = x2;
                    
                    x1 = x2;
                    
                    layer_area_available = 0;

                end
            end

        else
            layer_area_available = layer_area_available + bodyarea;
        end
        
        if (layersmade >= ylayers)
            % return if we've made enough layers, partly to sidestep
            % numerical issues
            cornernodes = [ lastinnernodes(2), botouternode, topouternode, lastinnernodes(1) ];
            return;
        end
        
        shoex = [ xcore + xcoil, shoex ];    
        shoey = [ ycoilshoe/2, shoey ];

        % gobble up the shoe area, creating layers as necessary
        if windingarea < (layer_area_available + shoearea)

            starttrapzind = 1;

            for ind = 2:numel (shoex)
                gobblearea = trapz (shoex(starttrapzind:ind), shoey(starttrapzind:ind)) + layer_area_available;
                if gobblearea >= windingarea
                    % create a layer link and coil label location

                    links = [ links; topshoenids(ind-1), botshoenids(ind-1) ];

                    coillabelloc = [ coillabelloc; lastlayerstartx + (shoex(ind) - lastlayerstartx)/2, 0 ];
                    layersmade = layersmade + 1;
                    lastlayerstartx = shoex(ind);
                    starttrapzind = ind;
                    layer_area_available = 0;
                end
            end

            if ylayers - layersmade == 1
                % put the remaining layer in the remaining area
                coillabelloc = [ coillabelloc; lastlayerstartx + (shoex(end) - lastlayerstartx)/2, 0 ];
            elseif layersmade ~= ylayers
                error ('layer construction error')
            end

        else
            if ylayers - layersmade == 1
                % put the remaining layer in the remaining area
                coillabelloc = [ coillabelloc; lastlayerstartx + (shoex(end) - lastlayerstartx)/2, 0 ];
            elseif layersmade ~= ylayers
                error ('layer construction error, ')
            end
        end

    
    end
    
    cornernodes = [ lastinnernodes(2), botouternode, topouternode, lastinnernodes(1) ];
                   
end


function [x, y, Qx, Qy] = shoecurvepoints (shoecontrolfrac, xcore, xcoil, xshoebase, xshoegap, ycoilbase, ycoilgap, yshoegap)
% creates a curve based on a quadratic Bezier curve with three control
% points
%
% 

    % get the slope at the start of the curve of the shoe
    cp1 = [xcore+xcoil, ycoilgap/2];
 
    % get the slope at the point where the curve ends at the shoe gap
    cp3 = [xcore+xcoil+xshoebase-xshoegap, yshoegap/2];
    
    % first get line equation
    leftslope = ((ycoilgap - ycoilbase)/2) / (xcoil);
    c = cp1(2) - leftslope * cp1(1);
    % get intercept
    yint = leftslope * cp3(1) + c;
    
    % calculate vector pointing from the right point to the intercept
    Vlr = [cp3(1), yint] - cp1;
    
    % control point lies at specified fraction along this vector
    cp2 = cp1 + shoecontrolfrac * Vlr;
    
    % construct the control points for the Bezier curve
    Px = [cp1(1), cp2(1), cp3(1)]; 
    Py = [cp1(2), cp2(2), cp3(2)];
    
    [Qx, Qy] = bezierpoints(Px,Py,100);
    
    % calculate the length of each curve segment
    [~,rho] = cart2pol (Qx(2:end) - Qx(1:end-1), Qy(2:end) - Qy(1:end-1));
    
    % calculate the total length of the curve
    curvelength = sum(rho);
    
    % choose appropriate number of points for femm to construct the curve,
    % either one every half mm, or at least 3 points
    npoints = max (3, min (ceil (curvelength/5e-4), 96));
    
    x = Qx(2:floor((numel(Qx)-2)/npoints):end-1);
    y = Qy(2:floor((numel(Qy)-2)/npoints):end-1);
    
end

function [x, y, Qx, Qy, m, c] = basecurvepoints (xcore, xcoilbase, xcoilbody, ycoilbase, ycoilgap)

    % get the start and end points of the curve, which are the first and
    % third control points
    cp1 = [xcore, 0];
    cp3 = [xcore+xcoilbase, ycoilbase/2];
    
    % first get line equation of slot side
    m = ((ycoilgap - ycoilbase)/2) / (xcoilbody);
    c = cp3(2) - m * cp3(1);
    % get intercept with y = xcore
    yint = m * xcore + c;
    
    % place the control point at the intercept. This should ensure that the
    % curve joins smoothly at both ends with the rest of the coil surface
    cp2 = [xcore, yint];
    
    % construct the control points for the Bezier curve
    Px = [cp1(1), cp2(1), cp3(1)]; 
    Py = [cp1(2), cp2(2), cp3(2)];
    
    [Qx, Qy] = bezierpoints(Px,Py,100);
    
    % calculate the length of each curve segment
    [~,rho] = cart2pol (Qx(2:end) - Qx(1:end-1), Qy(2:end) - Qy(1:end-1));
    
    % calculate the total length of the curve
    curvelength = sum(rho);
    
    % choose appropriate number of points for femm to construct the curve,
    % either one every half mm, or at least 3 points
    npoints = max (3, min (ceil (curvelength/5e-4), 96));
    
    x = Qx(2:floor((numel(Qx)-2)/npoints):end-1);
    y = Qy(2:floor((numel(Qy)-2)/npoints):end-1);

end


function y = mxplusc (m, c, x)

    y = m .* x + c;
    
end