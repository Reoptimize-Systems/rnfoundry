function [nodes, links, cornernodes, shoegaplabelloc, coillabelloc] ...
    = internalslotnodelinks(ycoil, yshoegap, xcore, xcoil, xshoebase, xshoegap, layers, tol)
% creates a set of node, links and label locations for the points making up
% the inside of a coil slot in an iron cored machine
%
% Syntax
%
% [nodes, links, cornernodes, shoelabelloc, shoegaplabelloc, coillabelloc] = ...
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
%

% Created by Richard Crozier 2012

    % there are no shoes so return an empty label
    shoegaplabelloc = [];
    
    if yshoegap > ycoil
        yshoegap = ycoil;
    end
        
    if xshoebase < tol
        % there is no shoe in this case, so just make a line for the
        % coil flush with the tooth surface
        
        nodes = [ xcore + xcoil, ycoil/2;
                  xcore + xcoil, -ycoil/2 ];
        
        links = [ 1, 0 ];

        topinnershoenode = 0;
        botinnershoenode = 1;
        
        topouternode = 0;
        botouternode = 1;
                
    else
        % there is a shoe on top of the tooth

        if yshoegap < tol
            % The shoe is joined in the middle seems odd, but we'll allow
            % for it here

            if xshoegap < tol
                
                % the tip of the shoe meets in a sharp point
                %
                % In this case three nodes are required to draw the shoe

                % the top shoe nodes
                nodes = [ xcore + xcoil, ycoil/2; ...
                          xcore + xcoil + xshoebase, ycoil/2; ...
                          xcore + xcoil + xshoebase, 0 ];
                      
                      
                links = [ 0, 2;
                          2, 1; ];
 
                % add the bottom shoe nodes
                nodes = [ nodes; ...
                          xcore + xcoil, -ycoil/2; ...
                          xcore + xcoil + xshoebase, -ycoil/2 ];
                      
                links = [ links;
                          2, 3;
                          4, 2 ];

                topinnershoenode = 0;
                botinnershoenode = 3; 
                
                topouternode = 1;
                botouternode = 4;

                if xshoebase > tol
                    % link the base of the shoes, the coil won't be in here
                    
                    links = [ links;
                              3, 0 ];
                          
                    shoegaplabelloc = tricenter( [xcore + xcoil + xshoebase, 0], ...
                                                 [xcore + xcoil, ycoil/2], ...
                                                 [xcore + xcoil, -ycoil/2], ...
                                                 'incenter', false)';
                    
                end

            else
                % there is a blunt edge on the shoe
                %
                % In this case 4 nodes are required to draw the shoe
                
                nodes = [ xcore + xcoil, ycoil/2; ...
                          xcore + xcoil + xshoebase, ycoil/2; ...
                          xcore + xcoil + xshoebase, 0; ...
                          xcore + xcoil + xshoegap, 0];
                      
                links = [ 2, 1;
                          2, 3;
                          3, 0];
                      
                % now draw the other shoe 
                nodes = [ nodes; 
                          nodes(2,1), -nodes(2,2); ...
                          nodes(1,1), -nodes(1,2); ];
                      
                links = [ links; ...
                          4, 2;
                          5, 3 ];
                      
                topinnershoenode = 0;
                botinnershoenode = 5; 
                
                topouternode = 1;
                botouternode = 4;
                
                if (xshoebase - xshoegap) > tol
                    % link the base of the shoes, the coil won't be in here
                    
                    links = [ links;
                              5, 0 ];
                          
                    shoegaplabelloc = tricenter( [xcore + xcoil + xshoebase - xshoegap, 0], ...
                                                 [xcore + xcoil, ycoil/2], ...
                                                 [xcore + xcoil, -ycoil/2], ...
                                                 'incenter', false)';
                    
                end
                
            end


        else
            % we must acount for the space between shoes in the model 
            if xshoegap < tol
                % the tip of the shoe ends in a sharp point

                % the top shoe nodes
                nodes = [ xcore + xcoil, ycoil/2; ...
                          xcore + xcoil + xshoebase, ycoil/2; ...
                          xcore + xcoil + xshoebase, yshoegap/2 ];
                      
                links = [ 1, 2;
                          2, 0 ];
                      
%                 shoelabelloc = tricenter(nodes(1,:), nodes(2,:), nodes(3,:), 'incenter', false)';    
                
                nodes = [ nodes; ...
                          xcore + xcoil, -ycoil/2; ...
                          xcore + xcoil + xshoebase, -ycoil/2; ...
                          xcore + xcoil + xshoebase, -yshoegap/2 ];
                      
                links = [ links; ...
                          4, 5;
                          3, 5 ];
                      
                links = [ links; 5, 2 ];

                topinnershoenode = 0;
                botinnershoenode = 3; 
                
                topouternode = 1;
                botouternode = 4;
                
                if xshoebase > tol
                    % link the base of the shoes, the coil won't be in here
                    
                    links = [ links;
                              3, 0 ];
                          
                    shoegaplabelloc = tricenter( [xcore + xcoil + xshoebase, 0], ...
                                                 [xcore + xcoil, ycoil/2], ...
                                                 [xcore + xcoil, -ycoil/2], ...
                                                 'incenter', false)';
                    
                end
                          
            else
                % there is a blunt edge on the shoe
                nodes = [ xcore + xcoil, ycoil/2; ...
                          xcore + xcoil + xshoebase, ycoil/2; ...
                          xcore + xcoil + xshoebase, yshoegap/2; ...
                          xcore + xcoil + xshoebase - xshoegap, yshoegap/2];
                      
                links = [ 2, 1;
                          2, 3;
                          0, 3];
                      
                % now draw the other shoe 
                nodes = [ nodes; 
                          nodes(:,1), -nodes(:,2); ];

                links = [ links; ...
                          5, 6;
                          6, 7;
                          7, 4 ];

                % make the shoe gap region
                links = [ links; 
                          6, 2; 
                          7, 3 ];
    
                % get the centre of the shoe gap
                shoegaplabelloc = rectcentre(nodes(8,:), nodes(3,:));
                       
                topinnershoenode = 0;
                botinnershoenode = 4; 
                
                topouternode = 1;
                botouternode = 5;
                
                if (xshoebase - xshoegap) > tol
                    % link the base of the shoes, the coil won't be in here
                    
                    links = [ links;
                              4, 0 ];
                          
                    shoegaplabelloc = [ shoegaplabelloc; ...
                                        tricenter( [xcore + xcoil + xshoebase - xshoegap, 0], ...
                                                   [xcore + xcoil, ycoil/2], ...
                                                   [xcore + xcoil, -ycoil/2], ...
                                                   'incenter', false)'; ...
                                      ];
                    
                end
            end


        end

    end

    
    % get the starting inner nodes for making links
    lastinnernodes = [ topinnershoenode, botinnershoenode ];
    
    % starting position for layers, the top of the outer layer
    layerpos = xcore + xcoil;
    
    for i = 1:layers
        % move the layerpos to the top of the next layer
        layerpos = layerpos - xcoil/layers;
        
        % get the number of the next node
        nextnode = size(nodes, 1);

        % create coil layer
        nodes = [ nodes; 
                  layerpos, ycoil/2; 
                  layerpos, -ycoil/2 ];

        links = [ links; ...
                  lastinnernodes(1),  nextnode; ...
                  lastinnernodes(2), nextnode+1; ...
                  nextnode + 1, nextnode ];
              
        % store the new coil part nodes 
        lastinnernodes = [ nextnode, nextnode + 1 ];      
        
        % get a suitible position for the coil label         
        coillabelloc(i,1:2) = rectcentre( [layerpos, -ycoil/2], ...
                                          [layerpos + xcoil/layers, ycoil/2]);
        
    end
    
    cornernodes = [ lastinnernodes(2), botouternode, topouternode, lastinnernodes(1) ];
                   
end