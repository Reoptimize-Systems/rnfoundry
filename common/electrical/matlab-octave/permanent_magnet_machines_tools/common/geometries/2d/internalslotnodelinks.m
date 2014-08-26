function [nodes, links, info] = internalslotnodelinks(ycoil, yshoegap, xcore, xcoil, xshoebase, xshoegap, ylayers, tol, varargin)
% creates a set of node, links and label locations for the points making up
% the inside of a coil slot in an iron cored machine
%
% Syntax
%
% [nodes, links, info] = ...
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
%  ycoil - If a single number, the internal slot width in y direction. If
%    a vector of two numbers, the first is the internal slot width at the
%    shoe end of the tooth, the second the internal slot width at the base
%    of the slot.
% 
%  yshoegap - distance between shoe tips in y direction. This can be any
%    dimension from zero up to ycoil at the shoe end.
% 
%  xcore - thickness of slot yoke in x direction (in practice this is the
%    distance of the base of the drawn slot shape from the y-axis) 
% 
%  xcoil - internal slot height in x direction, if layers are specified
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
%  In addition several options can be supplied as parameter-value pairs.
%  The possible optional arguments are:
%
%  CoilBaseFraction - fraction of the base at which curvature is present.
%    The slot base is curved with a bezier curve starting at xcore and
%    ending at an x position given by xcore + xcoil*CoilBaseFraction.
%    Defaults to 0.05 if not supplied.
%
%  ShoeCurveControlFrac - factor controlling the 'curvature' of the tooth
%    shoe, this is a value between 0 and 1. The exact effect of this number
%    is complex, and depends on the geometry of the slot. However, in
%    general a lower number results in a curve closer to a line draw
%    directly from the shoe base to the shoe gap, whicle higher numbers
%    aproximate a sharp right angle. Defaults to 0.5. 
%
%    N.B. the slot geometry affects this curve in the following way. If the
%    position of the shoe gap node is below the intercept of the line
%    formed by the edge of the slot and a vertical line at the shoe gap
%    node, the resulting curve will bend outward from the inside of the
%    slot. If the intercept is below the shoe gap node, the curve will bend
%    into the slot.
%
%  SplitX - if there are layers in the y direction, the slot can be split
%    into two in the x direction by setting this flag to true. Defaults to
%    false. If true coil label locations are provide from bottom to top.
%
% Output
%
% nodes - (n x 2) matrix of x and y coordinates of the nodes making up the
%  internal dimensions of the slots
%
% links - (n x 2) matrix of zero-based links between the nodes. The values
%   in links are the rows of 'nodes' counting from zero. 
% 
% info.cornernodes - vector containin the zero-based rows of the nodes at the
%   outer corners of the slot shape. This are supplied clockwise from the
%   bottom left, as shown in Figure 2, such that the vector contains
%
%   [ cn1, cn2, cn3, cn4 ]
%  
% info.shoegaplabelloc - (n x 2) matrix of label locations of any gaps created
%   by the shoe geometry. If yshoegap (and yshoebase) is greater than zero,
%   there will be an air-filled rectangular space between the shoe tips. A
%   suitable label location for this region is provided here. If the shoes
%   have a wider base than tip, there will be an additional trianglular or
%   trapezoidal gap between the top of the coil and the start of the space
%   between the shoe tips, and another suitable label location will be
%   provided. If xshoebase is zero (so there are no shoes), info.shoegaplabelloc
%   will be empty.
% 
% info.coillabelloc - (n x 2) matrix of coil label locations. The number of rows
%   will depend on the value of layers, with one label for every layer.
%
% info.vertlinkinds - vector of (one-based) indices of the links which should be 
%   converted to arcs when using internalslotnodelinks to make a radial 
%   design, e.g the outer tooth surface links.
%
% info.toothlinkinds - vector of indices of the links that make up the internal
%   surface of the tooth, not including the links creating divisions in the
%   cross-section of the slot of coils, or the marking the air region
%   between the shoes if there is a gap.
%

% Copyright Richard Crozier 2012-2014

    if nargin == 0 && nargout == 1
        % return function handles to subfunctions for testing purposes
        nodes = {@shoecurvepoints, @basecurvepoints, @inscurvepoints};
        return;
    end
    
    options.CoilBaseFraction = 0.05;
    options.ShoeCurveControlFrac = 0.5;
    options.SplitX = false;
    options.MaxBaseCurvePoints = 20;
    options.MaxShoeCurvePoints = 20;
    options.MinBaseCurvePoints = 5;
    options.MinShoeCurvePoints = 5;
    options.InsulationThickness = 0;
    options.YScale = 1;
    
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
    
    % initialise empty label locations
    info.shoegaplabelloc = [];
    info.inslabelloc = [];
    
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
    
    info.inslinkinds = [];
        
    if xshoebase < tol
        % there is no shoe in this case, so just make a line for the
        % coil flush with the tooth surface
        
        nodes = [ xcore + xcoilbase + xcoilbody, ycoilshoe/2;
                  xcore + xcoilbase + xcoilbody, -ycoilshoe/2;
                  xcore + xcoilbase + xcoilbody, 0 ];
        
        links = [ 1, 2;
                  2, 0 ];
                  
        nshoecurvepts = 0;
              
        shoex = [];
        shoey = [];
        topshoenids = [];
        botshoenids = [];

        topinnershoenode = 0;
        botinnershoenode = 1;
        
        topouternode = 0;
        botouternode = 1;
        
        midshoenode = 2;
        
        info.vertlinkinds = [1, 2];
        
        info.toothlinkinds = [];
                
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
                [shoex, shoey, shoeQx, shoeQy, shoePx, shoePy] ...
                               = shoecurvepoints ( options.ShoeCurveControlFrac, ...
                                                   [ options.MinShoeCurvePoints, options.MaxShoeCurvePoints ], ...
                                                   xcore, ...
                                                   xcoilbase, ...
                                                   xcoilbody, ...
                                                   xshoebase, ...
                                                   xshoegap, ...
                                                   ycoilbase, ...
                                                   ycoilshoe, ...
                                                   yshoegap );
                                 
                nshoecurvepts = numel(shoex);
                
                % the top shoe nodes
                nodes = [ xcore + xcoilbase + xcoilbody, ycoilshoe/2; ...
                          xcore + xcoilbase + xcoilbody + xshoebase, ycoilshoe/2; ...
                          xcore + xcoilbase + xcoilbody + xshoebase, 0; ...
                          xcore + xcoilbase + xcoilbody, -ycoilshoe/2; ...
                          xcore + xcoilbase + xcoilbody + xshoebase, -ycoilshoe/2; ...
                          shoex', shoey'; ...
                          shoex', -shoey' ];                       
                
                topshoenids = 5:4+nshoecurvepts;
                botshoenids = 5+nshoecurvepts:4+2*nshoecurvepts;
                
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
                
                info.vertlinkinds = [ 1, size(links,1)];

            else
                % there is a blunt edge on the shoe
                %
                % In this case 4 nodes are required to draw the shoe
                % corners
                
                % calculate intermediate points between the shoe and top of
                % coil region which create a curved surface
                [shoex, shoey, shoeQx, shoeQy, shoePx, shoePy] ...
                               = shoecurvepoints ( options.ShoeCurveControlFrac, ...
                                                   [ options.MinShoeCurvePoints, options.MaxShoeCurvePoints ], ...
                                                   xcore, ...
                                                   xcoilbase, ...
                                                   xcoilbody, ...
                                                   xshoebase, ...
                                                   xshoegap, ...
                                                   ycoilbase, ...
                                                   ycoilshoe, ...
                                                   yshoegap );

                nshoecurvepts = numel(shoex);
                
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

                topshoenids = 6:5+nshoecurvepts;
                botshoenids = 6+nshoecurvepts:5+2*nshoecurvepts;
                
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
                
                info.vertlinkinds = [ 1, size(links,1)];
                
                
            end
            
            % all the links are the tooth surface at this stage
            info.toothlinkinds = 1:size(links,1);

        else
            % we must acount for the space between shoes in the model 
            if xshoegap < tol
                % the tip of the shoe ends in a sharp point

                [shoex, shoey, shoeQx, shoeQy, shoePx, shoePy] ...
                               = shoecurvepoints ( options.ShoeCurveControlFrac, ...
                                                   [ options.MinShoeCurvePoints, options.MaxShoeCurvePoints ], ...
                                                   xcore, ...
                                                   xcoilbase, ...
                                                   xcoilbody, ...
                                                   xshoebase, ...
                                                   xshoegap, ...
                                                   ycoilbase, ...
                                                   ycoilshoe, ...
                                                   yshoegap );
                                            
                nshoecurvepts = numel(shoex);
                
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
                      
                topshoenids = 7:6+nshoecurvepts;
                botshoenids = 7+nshoecurvepts:6+2*nshoecurvepts;
                
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
                
                info.vertlinkinds = [ 1, size(links,1)-2];
                
                % tooth links are all links except the links closing the
                % air gap
                info.toothlinkinds = 1:size(links,1);
                info.toothlinkinds([size(links,1), size(links,1)-1]) = [];
                          
            else
                % blunt edged shoes with gap
                
                % calculate intermediate points between the shoe and top of
                % coil region which create a curved surface
                [shoex, shoey, shoeQx, shoeQy, shoePx, shoePy] ...
                               = shoecurvepoints ( options.ShoeCurveControlFrac, ...
                                                   [ options.MinShoeCurvePoints, options.MaxShoeCurvePoints ], ...
                                                   xcore, ...
                                                   xcoilbase, ...
                                                   xcoilbody, ...
                                                   xshoebase, ...
                                                   xshoegap, ...
                                                   ycoilbase, ...
                                                   ycoilshoe, ...
                                                   yshoegap );

                nshoecurvepts = numel(shoex);
                
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
                      
                topshoenids = 9:8+nshoecurvepts;
                botshoenids = 9+nshoecurvepts:8+2*nshoecurvepts;
                
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
                info.toothlinkinds = 1:size(links,1);
                
                % make the shoe gap region
                links = [ links; 
                          6, 2; 
                          7, 8;
                          8, 3 ];
    
                % get the centre of the shoe gap
                info.shoegaplabelloc = rectcentre(nodes(8,:), nodes(3,:));
                       
                topinnershoenode = 0;
                botinnershoenode = 4; 
                
                topouternode = 1;
                botouternode = 5;
                
                midshoenode = 8;
                
                info.vertlinkinds = [ 1, 2, size(links,1)-2, size(links,1)-1, size(links,1)-0 ];
                
            end

        end

    end
    
    % make the curve at the base of the slot
    [basex, basey, baseQx, baseQy, m, c, basePx, basePy] = basecurvepoints ( ...
                    [ options.MinBaseCurvePoints, options.MaxBaseCurvePoints ], ...
                    xcore, xcoilbase, xcoilbody, ycoilbase, ycoilshoe );
    
    nbasecurvepnts = numel (basex);
    
    % this will be the id of the next node we add, the (xcore, 0) node
    startbasenodeid = size (nodes, 1);
    
    nodes = [ nodes;
              xcore, 0; 
              basex', basey'; 
              basex', -basey';
              xcore+xcoilbase, ycoilbase/2;
              xcore+xcoilbase, -ycoilbase/2 ];
          
    endbasenodeid = [ size(nodes, 1) - 1, size(nodes, 1) - 2];

    linksize = size(links,1);
    
    links = [ links;
              startbasenodeid, startbasenodeid+1; 
              [(startbasenodeid+1:startbasenodeid+nbasecurvepnts-1)', (startbasenodeid+2:startbasenodeid+nbasecurvepnts)']
              startbasenodeid+nbasecurvepnts, endbasenodeid(2);
              startbasenodeid, startbasenodeid+nbasecurvepnts+1;
              [(startbasenodeid+nbasecurvepnts+1:startbasenodeid+2*nbasecurvepnts-1)', (startbasenodeid+nbasecurvepnts+2:startbasenodeid+2*nbasecurvepnts)'] 
              startbasenodeid+2*nbasecurvepnts, endbasenodeid(1)  ];
    
    % store the tooth surface link indices
    info.toothlinkinds = [info.toothlinkinds, linksize+1:size(links,1)];
    
    % get the starting inner nodes for making links
    lastinnernodes = [ topinnershoenode, botinnershoenode ];
    
    % starting position for layers, the top of the outer layer
    
    % calculate area of various coil regions, we will operate on one half
    % of the winding to simplify finding the area under curves etc
    basex = [xcore, basex, xcore+xcoilbase];
    basey = [0, basey, ycoilbase/2];
    topbasenids = [(startbasenodeid:startbasenodeid+nbasecurvepnts), endbasenodeid(2)];
	botbasenids = [startbasenodeid, (startbasenodeid+nbasecurvepnts+1:startbasenodeid+2*nbasecurvepnts), endbasenodeid(1)];
    
    if options.InsulationThickness > 0
    
        % reduce the intercept of the slot side curve, to match the insulation
        % thickness
        c = c - options.InsulationThickness/options.YScale;
        
        % get the number of existing links, useful for several calcs
        linksize = size (links,1);
        
        % create the shoe curve inulation points and links
    
        if xshoebase < tol
        
            nnodes = size (nodes,1);
            
            % add an inner insulation node at the 
            nodes = [ nodes; ...
                      xcore + xcoilbase + xcoilbody - options.InsulationThickness, mxplusc(m, c, xcore + xcoilbase + xcoilbody - options.InsulationThickness); ...
                      xcore + xcoilbase + xcoilbody - options.InsulationThickness, -mxplusc(m, c, xcore + xcoilbase + xcoilbody - options.InsulationThickness); ...
                      xcore + xcoilbase + xcoilbody - options.InsulationThickness, 0; ];
                          
            links = [ links;
                      [ 1, 2; 2, 0] + nnodes ];
                      
            nshoecurvepts = 0;
                  
            shoex = [];
            shoey = [];
            topshoenids = [];
            botshoenids = [];

            topinnershoenode = 0 + nnodes;
            botinnershoenode = 1 + nnodes;
            
            topouternode = 0 + nnodes;
            botouternode = 1 + nnodes;
            
            midshoenode = 2 + nnodes;
            
            info.vertlinkinds = [info.vertlinkinds, [1, 2] + linksize];
        
        else
        
            % get the shoe insulation curve points
            [shoex, shoey, shoeQx, shoeQy] = inscurvepoints ( ...
                            [ options.MinShoeCurvePoints, options.MaxShoeCurvePoints ], ...
                            shoeQx, shoeQy, shoePx, shoePy, ...
                            options.InsulationThickness, options.YScale );
                  
            
%            shoex = shoex(shoey(1)*options.YScale - shoey*options.YScale >= options.InsulationThickness);
%            shoey = shoey(shoey(1)*options.YScale - shoey*options.YScale >= options.InsulationThickness);
            
            crossinginds = find (shoey < -shoey);
            
            if ~isempty (crossinginds)
                
                shoex = shoex (1:crossinginds(1)-1);
                shoey = shoey (1:crossinginds(1)-1);
                
                nodes = [ nodes; ...
                          shoex(end), 0 ; ];
            else
                nodes = [ nodes; ...
                          xcore + xcoilbase + xcoilbody + xshoebase - xshoegap - options.InsulationThickness, 0 ; ];
            end
            
            if isempty (shoex)
                shoex = xcore + xcoilbase + xcoilbody - options.InsulationThickness;
                shoey = mxplusc(m, c, xcore + xcoilbase + xcoilbody - options.InsulationThickness);
            end
                            
            nshoecurvepts = numel (shoex);
                      
            nnodes = size(nodes, 1);
            
            nodes = [ nodes; ...
                      shoex', shoey';
                      shoex', -shoey' ];
                      
            topshoenids = nnodes:nnodes-1+nshoecurvepts;
            botshoenids = topshoenids(end)+1:topshoenids(end)+nshoecurvepts;
            midshoenode = nnodes - 1;
            topinnershoenode = topshoenids(1);
            botinnershoenode = botshoenids(1);
            
            links = [ links;
                      [topshoenids(1:end-1)', topshoenids(2:end)' ]; ... % join the curve nodes
                      topshoenids(end), midshoenode; ...
                      midshoenode, botshoenids(end); ...
                      [botshoenids(1:end-1)', botshoenids(2:end)' ]; ];
                  
            info.vertlinkinds = [info.vertlinkinds, size(links,1)-numel(botshoenids), size(links,1)-numel(botshoenids)+1];
            
        end
        
        % add the new insulation links to the list
        info.inslinkinds = [ info.inslinkinds, linksize+1:size(links,1) ];
        
        linksize = size (links,1);
        
        % we must link up the tooth sides and replace basex and basey with the 
        % internal insulation basex and basey 
        links = [ links; 
                  lastinnernodes(1), endbasenodeid(2);
                  lastinnernodes(2),  endbasenodeid(1) ];
                  
        info.toothlinkinds = [info.toothlinkinds, linksize+1:size(links,1)];
        
        % get the insulation curve points
        [basex, basey, baseQx, baseQy] = inscurvepoints ( ...
                        [ options.MinBaseCurvePoints, options.MaxBaseCurvePoints ], ...
                        baseQx, baseQy, basePx, basePy, ...
                        options.InsulationThickness, options.YScale );
                        
        basex = basex(basey(end)*options.YScale - basey*options.YScale >= options.InsulationThickness);
        basey = basey(basey(end)*options.YScale - basey*options.YScale >= options.InsulationThickness);
                        
        crossinginds = find (basey > mxplusc(m, c, basex));
        
        if ~isempty (crossinginds)
            basex = basex (1:crossinginds(1)-1);
            basey = basey (1:crossinginds(1)-1);
        end
        
        crossinginds = find (basey <= 0);
        
        if ~isempty (crossinginds)
            basex = basex (crossinginds(end)+1:end);
            basey = basey (crossinginds(end)+1:end);
        end
            
        nbasecurvepnts = numel (basex);
        
        % this will be the id of the next node we add, the
        %  (xcore + options.InsulationThickness, 0) node
        startbasenodeid = size (nodes, 1);
        
        nodes = [ nodes;
                  basex(1), 0; 
                  basex', basey'; 
                  basex', -basey'; ];
              
        endbasenodeid = [ startbasenodeid+2*nbasecurvepnts, startbasenodeid+nbasecurvepnts ];

        linksize = size(links,1);
        
        % add the insulation base links
        links = [ links;
                  startbasenodeid, startbasenodeid+1; 
                  [ (startbasenodeid+1:startbasenodeid+nbasecurvepnts-1)', ...
                    (startbasenodeid+2:startbasenodeid+nbasecurvepnts)']
                  startbasenodeid, startbasenodeid+nbasecurvepnts+1;
                  [ (startbasenodeid+nbasecurvepnts+1:startbasenodeid+2*nbasecurvepnts-1)', ...
                    (startbasenodeid+nbasecurvepnts+2:startbasenodeid+2*nbasecurvepnts)'] ...
                ];
                
        info.inslinkinds = [ info.inslinkinds, linksize+1:size(links,1) ];
                
        basex = [xcore + options.InsulationThickness, basex];
        basey = [0, basey];
        topbasenids = (startbasenodeid:startbasenodeid+nbasecurvepnts);
        botbasenids = [startbasenodeid, (startbasenodeid+nbasecurvepnts+1:startbasenodeid+2*nbasecurvepnts)];
        
        info.inslabelloc = [xcore + options.InsulationThickness/2, 0];
        
        lastinnernodes = [topinnershoenode, botinnershoenode];
        
    end
    
    basearea = trapz (basex, basey);
    bodyarea = trapzarea ( ycoilbase, ycoilshoe, xcoilbody ) / 2;
    shoearea = trapz ([ xcore+xcoil, shoex, xcore+xcoil+xshoebase-xshoegap ], ...
                      [ ycoilbase/2, shoey, yshoegap ] );
    
    info.totalarea = basearea + bodyarea + shoearea;
    
    % calculate the area of one winding layer, given the number of layers
    info.windingarea = info.totalarea / ylayers;
    
    linksize = size(links,1);
    
    if ylayers == 1
        % simple case, there is only one layer per slot, winding fills the
        % entire slot, but with the option of splitting into two layers in
        % the alternate direction so coils lie side-by-side one another in
        % the slot instead of on top of one another
        
        if options.SplitX
            
            % add a link creating the sides and a horizontal split
            links = [ links; 
                      startbasenodeid, midshoenode;
                      lastinnernodes(1), endbasenodeid(2);
                      lastinnernodes(2),  endbasenodeid(1) ];
            
            % get a suitible position for the coil labels
            info.coillabelloc = [ xcore + xcoilbase + xcoilbody/2, ycoilshoe/4;
                             xcore + xcoilbase + xcoilbody/2, -ycoilshoe/4 ];
        else
            % just add links making the sides
            links = [ links; 
                      lastinnernodes(1), endbasenodeid(2);
                      lastinnernodes(2),  endbasenodeid(1) ];
                  
            % get a suitible position for the coil label         
            info.coillabelloc = [ xcore + (xcoilbase + xcoilbody)/2, 0 ];
        end
        
        % update the link lists as appropriate
        if options.InsulationThickness > 0
            info.inslinkinds = [info.inslinkinds, size(links, 1) - 1, size(links, 1)];
        else
            info.toothlinkinds = [info.toothlinkinds, size(links, 1) - 1, size(links, 1)];
        end
        
        info.cornernodes = [ lastinnernodes(2), botouternode, topouternode, lastinnernodes(1) ];
        
    else
    
        if options.SplitX
            error ('RENEWNET:internalslotnodelinks:badoption', ...
                   'Option ''SplitX'' is incompatible with coillayers > 1');
        end
        
        layersmade = 0;
        
        info.coillabelloc = [];
        
        basemidnode = startbasenodeid;
        
        layer_area_available = 0;
        
        lastlayerstartx = xcore;
        
        if info.windingarea < basearea
            % gobble up the base area, creating layers as necessary
            starttrapzind = 1;
                    
            for ind = 2:numel (basex)
                gobblearea = trapz (basex(starttrapzind:ind), basey(starttrapzind:ind));
                if gobblearea >= info.windingarea
                    % create a layer link and coil label location
                    
                    links = [ links; topbasenids(ind), botbasenids(ind) ];
                    
                    info.vertlinkinds = [info.vertlinkinds, size(links, 1)];
                    
                    info.coillabelloc = [ info.coillabelloc; lastlayerstartx + (basex(ind) - lastlayerstartx)/2, 0 ];
                    layersmade = layersmade + 1;
                    
                    layer_area_available = info.windingarea - gobblearea;
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
            info.cornernodes = [ lastinnernodes(2), botouternode, topouternode, lastinnernodes(1) ];
            return;
        end

        if info.windingarea < (layer_area_available + bodyarea)
            % we calculate the required position for x by finding the roots
            % of a quadratic obtained from the integral of the staight line
            % making up the straight side of the coil
            x1 = xcore + xcoilbase;
            
            lastlayernodeids = [topbasenids(end), botbasenids(end)];
            
            while 1
                
                bodyarealeft = trapzarea ( mxplusc (m, c, x1)*2, ...
                                           mxplusc (m, c, xcore + xcoilbase + xcoilbody)*2, ...
                                           (xcore + xcoil) - x1 ) / 2;

                bodyarealeft = bodyarealeft + layer_area_available;
                
                if (info.windingarea > bodyarealeft) || (layersmade >= ylayers) ...

                    % link up the sides
                    links = [ links;
                              lastlayernodeids(1), topinnershoenode;
                              lastlayernodeids(2), botinnershoenode; ];
                    
                    % update the link lists
                    if options.InsulationThickness > 0
                        info.inslinkinds = [info.inslinkinds, size(links, 1) - 1, size(links, 1)];
                    else
                        info.toothlinkinds = [info.toothlinkinds, size(links, 1) - 1, size(links, 1)];
                    end
                    
                    % store the remaining area in the body available to
                    % a layer
                    layer_area_available = bodyarealeft;
                    
                    % break out of the while loop
                    break;
                    
                elseif ((xshoebase < tol) && (layersmade == (ylayers - 1))) ...
                        || ((xshoebase - xshoegap) < tol && (layersmade == (ylayers - 1)))
                    
                    % link up the sides
                    links = [ links;
                              lastlayernodeids(1), topinnershoenode;
                              lastlayernodeids(2), botinnershoenode; ];

                    % update the link lists
                    if options.InsulationThickness > 0
                        info.inslinkinds = [info.inslinkinds, size(links, 1) - 1, size(links, 1)];
                    else
                        info.toothlinkinds = [info.toothlinkinds, size(links, 1) - 1, size(links, 1)];
                    end

                    % put the remaining layer in the space between the last area
                    % boundary and the top
                    info.coillabelloc = [ info.coillabelloc; lastlayerstartx + ((xcore+xcoilbody+xcoilbase) - lastlayerstartx)/2, 0 ];

                    % return as there is no shoe area to gobble
                    info.cornernodes = [ lastinnernodes(2), botouternode, topouternode, lastinnernodes(1) ];
                    return;
                    
                else
                    % calculate the area of the body required to make up the
                    % remaining required winding area
                    At = info.windingarea - layer_area_available;

                    y1 = mxplusc (m, c, x1);
                    
                    % find the roots of the polynomial to get the position of x
                    arearoots = roots ([ m/2, ...
                                         (y1 + c -m*x1)/2, ...
                                         (-y1*x1 - c*x1)/2 - At ]);
                                     
                    % there will be two valid roots, choose the one which
                    % gives positive y2
                    y2 = mxplusc(m,c,arearoots);

                    x2 = arearoots(y2>0);

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
                          
                    info.vertlinkinds = [info.vertlinkinds, size(links, 1)];
                          
                    % update the link lists
                    if options.InsulationThickness > 0
                        info.inslinkinds = [info.inslinkinds, size(links, 1) - 1, size(links, 1)];
                    else
                        info.toothlinkinds = [info.toothlinkinds, size(links, 1) - 1, size(links, 1)];
                    end

                    info.coillabelloc = [ info.coillabelloc; lastlayerstartx + (x2-lastlayerstartx)/2, 0 ];
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
            info.cornernodes = [ lastinnernodes(2), botouternode, topouternode, lastinnernodes(1) ];
            return;
        end
        
        shoex = [ xcore + xcoil, shoex ];
        shoey = [ ycoilshoe/2 - options.InsulationThickness, shoey ];


        % gobble up the shoe area, creating layers as necessary
        if info.windingarea < (layer_area_available + shoearea)

            starttrapzind = 1;

            for ind = 2:numel (shoex)
                gobblearea = trapz (shoex(starttrapzind:ind), shoey(starttrapzind:ind)) + layer_area_available;
                if gobblearea >= info.windingarea
                    % create a layer link and coil label location

                    links = [ links; topshoenids(ind-1), botshoenids(ind-1) ];
                    
                    info.vertlinkinds = [info.vertlinkinds, size(links, 1)];

                    info.coillabelloc = [ info.coillabelloc; lastlayerstartx + (shoex(ind) - lastlayerstartx)/2, 0 ];
                    layersmade = layersmade + 1;
                    lastlayerstartx = shoex(ind);
                    starttrapzind = ind;
                    layer_area_available = 0;
                end
            end

            if ylayers - layersmade == 1
                % put the remaining layer in the remaining area
                info.coillabelloc = [ info.coillabelloc; lastlayerstartx + (shoex(end) - lastlayerstartx)/2, 0 ];
            elseif layersmade ~= ylayers
                error ('layer construction error')
            end

        else
            if ylayers - layersmade == 1
                % put the remaining layer in the remaining area
                info.coillabelloc = [ info.coillabelloc; lastlayerstartx + (shoex(end) - lastlayerstartx)/2, 0 ];
            elseif layersmade ~= ylayers
                error ('layer construction error, ')
            end
        end

    
    end
    
    info.cornernodes = [ lastinnernodes(2), botouternode, topouternode, lastinnernodes(1) ];
                   
end


function [x, y, Qx, Qy, Px, Py] = shoecurvepoints (shoecontrolfrac, minmaxpnts, xcore, xcoilbase, xcoilbody, xshoebase, xshoegap, ycoilbase, ycoilgap, yshoegap, tinsulate)
% creates a curve based on a quadratic Bezier curve with three control
% points
%
% 

    % get the slope at the start of the curve of the shoe
    cp1 = [xcore+xcoilbase+xcoilbody, ycoilgap/2];
 
    % get the slope at the point where the curve ends at the shoe gap
    cp3 = [xcore+xcoilbase+xcoilbody+xshoebase-xshoegap, yshoegap/2];
    
    % first get line equation
    leftslope = ((ycoilgap - ycoilbase)/2) / (xcoilbody);
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
    
    % calculate the total length of the curve
    curvelength = sum(magn ([Qx(2:end); Qy(2:end)] - [Qx(1:end-1); Qy(1:end-1)]));
    
    % choose appropriate number of points for femm to construct the curve,
    % either one every half mm, or at least 3 points
    npoints = max (minmaxpnts(1), min (ceil (curvelength/5e-4), minmaxpnts(2)));
    
    x = Qx(2:floor((numel(Qx)-2)/npoints):end-1);
    y = Qy(2:floor((numel(Qy)-2)/npoints):end-1);
    
end

function [insx, insy, insQx, insQy] = inscurvepoints (minmaxpnts, Qx, Qy, Px, Py, tins, yscale)
% creates a curve slightly inset from an existing bezier curve
 
    % calculate the total length of the curve
    curvelength = [0,0];
    
    sectionlengths = magn ([Qx(2:end); Qy(2:end)] - [Qx(1:end-1); Qy(1:end-1)]);
    
    startinsind = 2;
    endinsind = numel (sectionlengths);
    for ind = 1:numel (sectionlengths)
    
        curvelength = [ curvelength(1) + sectionlengths(ind), ...
                        curvelength(2) + sectionlengths(end-ind+1) ];
        
        if curvelength(1) <= tins
            startinsind = ind+1;
        end
        
        if curvelength(2) <= tins
            endinsind = numel(Qx) - ind;
        end
    
    end
    
    curvelength = sum (sectionlengths);
    
    % choose appropriate number of points for femm to construct the curve,
    % either one every half mm, or at least 3 points
    npoints = max (minmaxpnts(1), min (ceil (curvelength/5e-4), minmaxpnts(2)));
    
    % create a rotation matrix, to rotate the tangents 90 degrees
    rotangle = -tau/4;

    rotM = [ cos(rotangle)  sin(rotangle); 
             -sin(rotangle) cos(rotangle) ];
      
    n = numel (Qx);
    dt = 1/n;
    t = (0:n-1) * dt;

    P = [ Px; Py ];

    % the derivative of a qudratic bezier curve
    dQ = bsxfun (@times, 2 .* (1 - t), (P(:,2) - P(:,1))) ...
         + bsxfun (@times, 2 .* t, (P(:,3) - P(:,2)));

    normdQ = unit([ dQ(2,:); -dQ(1,:)]);
    
    normalshift = normdQ * tins;    
             
    % get locations at Qx and Qy, but shifted inward along a line perpendicular
    % to the tangent by the insulation thickness
    insQx = Qx + normalshift(1,:);
    insQy = Qy + normalshift(2,:)/yscale;
    
    insx = insQx(startinsind:floor((numel(insQx)-2)/npoints):endinsind);
    insy = insQy(startinsind:floor((numel(insQy)-2)/npoints):endinsind);
    
end

function [x, y, Qx, Qy, m, c, Px, Py] = basecurvepoints (minmaxpnts, xcore, xcoilbase, xcoilbody, ycoilbase, ycoilgap)

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
    cp2 = [xcore, max(ycoilbase/100,yint)];
    
    % construct the control points for the Bezier curve
    Px = [cp1(1), cp2(1), cp3(1)]; 
    Py = [cp1(2), cp2(2), cp3(2)];
    
    [Qx, Qy] = bezierpoints(Px,Py,100);
    
    % calculate the total length of the curve
    curvelength = sum(magn ([Qx(2:end); Qy(2:end)] - [Qx(1:end-1); Qy(1:end-1)]));
    
    % choose appropriate number of points for femm to construct the curve,
    % either one every half mm, or at least 5 points
    npoints = max (minmaxpnts(1), min (ceil (curvelength/5e-4), minmaxpnts(2)));
    
    x = Qx(2:floor((numel(Qx)-2)/npoints):end-1);
    y = Qy(2:floor((numel(Qy)-2)/npoints):end-1);

end


function y = mxplusc (m, c, x)

    y = m .* x + c;
    
end