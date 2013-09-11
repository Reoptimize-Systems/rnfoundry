function [nodes, links, cornernodes, shoelabelloc, shoegaplabelloc, coillabelloc] = internalslotnodelinks(ycoil, yshoegap, xcore, xcoil, xshoebase, xshoegap, layers, tol)
% creates a set of node, links and label locations for the points making up
% the inside of a coil slot in an iron cored machine
%
% Syntax
%
% [nodes, links, cornernodes, shoelabelloc, shoegaplabelloc, coillabelloc] = ...
%        internalslotnodelinks(ycoil, yshoegap, xcore, xcoil, xshoebase, xshoegap, layers, tol)
%
% Input
%
%  ycoil, 
% 
%  yshoegap, 
% 
%  xcore, 
% 
%  xcoil, 
% 
%  xshoebase, 
% 
%  xshoegap, 
% 
%  layers, 
% 
%  tol
%

    % there are no shoes so return an empty label
    shoelabelloc = [];
    shoegaplabelloc = [];
        
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

%                 shoelabelloc = tricenter(nodes(1,:), nodes(2,:), nodes(3,:), 'incenter', false)'; 
                
                % add the bottom shoe nodes
                nodes = [ nodes; ...
                          xcore + xcoil, -ycoil/2; ...
                          xcore + xcoil + xshoebase, -ycoil/2 ];
                      
                links = [ links;
                          2, 3;
                          4, 2 ];
                      
%                 shoelabelloc(2,:) = [shoelabelloc(1,1), -shoelabelloc(1,2)];

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
                      
                
%                 geom = polygeom( nodes(:,1), nodes(:,2) ) ; 
%                 shoelabelloc = geom(2:3);
                
                % now draw the other shoe 
                nodes = [ nodes; 
                          nodes(2,1), -nodes(2,2); ...
                          nodes(1,1), -nodes(1,2); ];
                      
                links = [ links; ...
                          4, 2;
                          5, 3 ];
                      
%                 shoelabelloc(2,:) = [shoelabelloc(1,1), -shoelabelloc(1,2)];
   
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
                      
%                 shoelabelloc(2,:) = [shoelabelloc(1,1), -shoelabelloc(1,2)];

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
                      
                
%                 geom = polygeom( nodes(:,1), nodes(:,2) ) ; 
%                 shoelabelloc = geom(2:3);
                
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
                      
%                 shoelabelloc(2,:) = [shoelabelloc(1,1), -shoelabelloc(1,2)];
                      
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