function [FemmProblem, wrapperthickness, innercentres, outercentres, nodeids, linktb] = ...
    wrappedannularsecmagaperiodic(FemmProblem, thetapole, thetamag, rmag, roffset, pos, wrapperthickness, varargin)
% creates a FemmProblem geometry of a radial section containing two magnets
% optionally wrapped with additional layers and with periodic edges
%
%
% Syntax
%
% [FemmProblem, wrapperthickness, innercentres, outercentres, nodeids] = ...
%   wrappedannularsecmagaperiodic(FemmProblem, thetapole, thetamag, rmag, ...
%                                 roffset, pos, wrapperthickness)
%
% [...] = wrappedannularsecmagaperiodic(..., 'param', value, ...)
%
% Description
%
% wrappedannularsecmagaperiodic creates a periodic geometry of two magnets
% with spaces in between, with a base position shown in the figure below.
% In addition, any number of annular sectors can be added either to inside
% or outside of the main region (like wrappers for the main region).
%                       
%             ********
%      *******        *                                    
%    **   *             *                                              
%     *     *             * 
%       *    ***************                                            
%        *    *             *                                              
%          *    *             *                                                                                          
%           *    *             *                                              
%             *    *    Mag 2    *                                                                                     
%              *    *             *                                              
%               *    *             *                                              
%                *    *             *            
%                 *    ***************                                                                                          
%                  /    *             *                                              
%                 / *    *             *  .........................
%                /   *    *             *                ^
%               /    *    *             *                 :
%              /      *    *************** ..^.......      :                            
%             /       *    *             *   :             :
%            /         *    *             *   :             :                                                                          
%           /          *    *             *   : thetamag    :                                    
%          /            *   *    Mag 1    *    :             : thetapole             
%   single internal     *    *             *   :             :                                      
%   wrapper example     *    *             *   :             :            
%                       *    *             * . v.........     :                                       
%                        *    ***************                 :                                       
%                        *    *             *                 :               
%                        *    *             *                 v                
%  x                     ******************** ..............................                                       
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
%  wrapperthickness - either an (n x 2) matrix or column vector of wrapper
%    thicknesses. If an (n x 2) matrix. the first column specifies the
%    thickness of any desired wrapper on the left of the magnets, and the
%    second column the thickness of wrappers on the right hand side. The
%    wrappers are added moving progressively further from the magnet
%    position (either to the left or right) down the rows of the matrix.
%    Wrappers with thicknesses less than a tolerance are not added. The
%    default tolerance is 1e-5, but this value can be changed using the
%    appropriate optional parameter value pair (see below). If
%    wrapperthickness is a column vector the same thicknesses are used on
%    both sides.
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
    Inputs.WrapperGroup = 0;
    Inputs.Tol = 1e-5;
    Inputs.MeshSize = -1;

    % parse the input arguments
    Inputs = parse_pv_pairs(Inputs, varargin);
    
    if isnumeric(Inputs.MagDirections) && isvector(Inputs.MagDirections)
        Inputs.MagDirections = {Inputs.MagDirections(1), Inputs.MagDirections(2)};
    end
    
    if (abs(thetapole - pi)*roffset) <= Inputs.Tol
        linktb = true;
    else
        linktb = false;
    end
    
    % remove input fields specific to this function and pass the
    % remaining fields to planarrectmagaperiodic as a series of p-v pairs
    planarrectmagaperiodicInputs = rmfield(Inputs, {'WrapperGroup'});
    
    % Convert the inputs structure to a series of p-v pairs
    planarrectmagaperiodicInputs = struct2pvpairs(planarrectmagaperiodicInputs);
    
    % first draw periodic magnet regions
    [FemmProblem, nodes, nodeids, links, magblockinds] = ...
        annularsecmagaperiodic(FemmProblem, thetapole, thetamag, rmag, roffset, pos, planarrectmagaperiodicInputs{:});
    
    if linktb
        % remove the segment linking the nodes at the 'top' of the drawing
        FemmProblem.Segments(numel(FemmProblem.Segments)) = [];
        % change the segments linking to the last two nodes to be added to
        % be linked instead to the first two nodes
        seglinks = getseglinks_mfemm(FemmProblem);
        for ind = 1:size(seglinks,1)
            if seglinks(ind,1) == nodeids(end-1)
                FemmProblem.Segments(ind).n0 = nodeids(1);
            elseif seglinks(ind,2) == nodeids(end-1)
                FemmProblem.Segments(ind).n1 = nodeids(1);
            elseif seglinks(ind,1) == nodeids(end)
                FemmProblem.Segments(ind).n0 = nodeids(2);
            elseif seglinks(ind,2) == nodeids(end)
                FemmProblem.Segments(ind).n1 = nodeids(2);
            end
        end
        seglinks = getarclinks_mfemm(FemmProblem);
        for ind = 1:size(seglinks,1)
            if seglinks(ind,1) == nodeids(end-1)
                FemmProblem.ArcSegments(ind).n0 = nodeids(1);
            elseif seglinks(ind,2) == nodeids(end-1)
                FemmProblem.ArcSegments(ind).n1 = nodeids(1);
            elseif seglinks(ind,1) == nodeids(end)
                FemmProblem.ArcSegments(ind).n0 = nodeids(2);
            elseif seglinks(ind,2) == nodeids(end)
                FemmProblem.ArcSegments(ind).n1 = nodeids(2);
            end
        end
        % remove the unnecessary final 2 nodes
        FemmProblem.Nodes(end-1:end) = [];
    end
    
    % Now correct the magnet angles which are specified as numeric values,
    % and therefore represent angles relative to a normal pointing in the
    % direction of the magnet block label
    for ind = 1:size(magblockinds, 1)
        if ~isempty(FemmProblem.BlockLabels(magblockinds(ind)).MagDir) ...
                && isscalar(FemmProblem.BlockLabels(magblockinds(ind)).MagDir) ...
                && isempty(FemmProblem.BlockLabels(magblockinds(ind)).MagDirFctn)
            
            [maglabeltheta, ~] = cart2pol(FemmProblem.BlockLabels(magblockinds(ind)).Coords(1), ...
                                    FemmProblem.BlockLabels(magblockinds(ind)).Coords(2));
                
            FemmProblem.BlockLabels(magblockinds(ind)).MagDir = ...
                FemmProblem.BlockLabels(magblockinds(ind)).MagDir + rad2deg(maglabeltheta);
        
        end
    end
    
    % now add the back iron nodes and links, depending on their thicknesses
    if size(wrapperthickness,2) < 2
        wrapperthickness = [wrapperthickness, wrapperthickness];
    elseif size(wrapperthickness,2) > 2
        error('wrapperthickness must be a scaler or a (1 x 2) vector or (n x 2) matrix')
    end
    
%     if anyalldims(wrapperthickness < 0)
%         error('wrapper thicknesses must all be greater that 0')
%     end
    
    elcount = elementcount_mfemm(FemmProblem);
    
    innercentres = zeros(size(wrapperthickness)) * NaN;
    
    % wrapperthickness(1,1) is the first inner region thickness
    if wrapperthickness(1,1) > Inputs.Tol
        % add the nodes and segments for the inner side
        %
        %                        |   /
        %  wrapperthickness      |  /
        % * sin(2*thetapole)     | /. 2*thetapole
        %                        |/__.__
        %
        %                     wrapperthickness 
        %                     * cos(2*thetapole)
        %
        
        if wrapperthickness(1,1) > roffset
            error('wrapper thickness cannot be greater than magnet inner radius.');
        end
        
        innerrad = roffset - rmag/2 - wrapperthickness(1,1);
        
        % First node is to left of first node in 'nodes' matrix. this is at
        % the bottom of the sim
        [FemmProblem, nodeinds, botnodeid] = addnodes_mfemm(FemmProblem, ...
                            innerrad, ...
                            0, ...
                            'InGroup', Inputs.WrapperGroup);
        
        
        
        if linktb
            lastbotnodeid = elcount.NNodes - size(nodes,1) + 2;
            topnodeid = botnodeid;
            botboundarymarker = '';
        else
            
            % Second node is to left of penultimate node in 'nodes' matrix
            [x,y] = pol2cart(2*thetapole, innerrad);
            [FemmProblem, nodeinds, topnodeid] = addnodes_mfemm(FemmProblem, ...
                                 x, ...
                                 y, ...
                                'InGroup', Inputs.WrapperGroup);

            % add a new periodic boundary for the top and bottom of the region
            [FemmProblem, boundind] = addboundaryprop_mfemm(FemmProblem, 'Left Wrap Annular Sec Mags Periodic', 4);

            botboundarymarker = FemmProblem.BoundaryProps(boundind).Name;
            
            % Seg with Periodic boundary at top
            FemmProblem = addsegments_mfemm(FemmProblem, ...
                                            elcount.NNodes - 2, ...
                                            topnodeid, ...
                                            'BoundaryMarker', FemmProblem.BoundaryProps(boundind).Name, ...
                                            'InGroup', Inputs.WrapperGroup);
                                        
            lastbotnodeid = elcount.NNodes - size(nodes,1);
            
        end
        
        % Seg at bottom
        FemmProblem = addsegments_mfemm(FemmProblem, ...
                                        lastbotnodeid, ...
                                        botnodeid, ...
                                        'BoundaryMarker', botboundarymarker, ...
                                        'InGroup', Inputs.WrapperGroup);


        % Add a node at the mid-point of the wrapper, the purpose of this
        % is to ensure an arc segmetn is never more than 180 degrees which
        % causes problems
        [x,y] = pol2cart(thetapole, innerrad);
        [FemmProblem, nodeinds, midnodeid] = addnodes_mfemm(FemmProblem, ...
                                        x, ...
                                        y, ...
                                        'InGroup', Inputs.WrapperGroup);
                        
        % Seg joining top and bottom (the two most recently added nodes)
        FemmProblem = addarcsegments_mfemm( FemmProblem, ...
                                            botnodeid, ...
                                            midnodeid, ...
                                            rad2deg(thetapole), ...
                                            'InGroup', Inputs.WrapperGroup );
                                        
        % Seg joining top and bottom (the two most recently added nodes)
        FemmProblem = addarcsegments_mfemm( FemmProblem, ...
                                            midnodeid, ...
                                            topnodeid, ...
                                            rad2deg(thetapole), ...
                                            'InGroup', Inputs.WrapperGroup );
                    
        [innercentres(1,1), innercentres(1,2)] = pol2cart(thetapole, innerrad+wrapperthickness(1,1)/2);
        
    else
        % Set the region thickness to be exactly zero so this can be tested
        % later
        wrapperthickness(1,1) = 0;
    end
    
    % now add all subsequent inner wrappers
    for i = 2:size(wrapperthickness, 1)
        
        if wrapperthickness(i,1) > Inputs.Tol
            
            innerrad = innerrad - wrapperthickness(i,1);
        
            if innerrad < Inputs.Tol
                error('Inner wrapper radii must all be greater than tolerance.')
            end
            
            lastbotnodeid = botnodeid;
            lasttopnodeid = topnodeid;
            
            % First node is to left of first node in 'nodes' matrix. this is at
            % the bottom of the sim
            [FemmProblem, nodeinds, botnodeid] = addnodes_mfemm(FemmProblem, ...
                                innerrad, ...
                                0, ...
                                'InGroup', Inputs.WrapperGroup);
   
            if linktb
                topnodeid = botnodeid;
                botboundarymarker = '';
            else
                
                % Second node is to left of penultimate node in 'nodes' matrix
                [x,y] = pol2cart(2*thetapole, innerrad);
                [FemmProblem, nodeinds, topnodeid] = addnodes_mfemm(FemmProblem, ...
                                     x, ...
                                     y, ...
                                    'InGroup', Inputs.WrapperGroup);

                % add a new periodic boundary for the top and bottom of the region
                [FemmProblem, boundind] = addboundaryprop_mfemm(FemmProblem, 'Left Wrap Annular Sec Mags Periodic', 4);

                botboundarymarker = FemmProblem.BoundaryProps(boundind).Name;
                
                % Seg with Periodic boundary at top
                FemmProblem = addsegments_mfemm(FemmProblem, ...
                                                lasttopnodeid, ...
                                                topnodeid, ...
                                                'BoundaryMarker', botboundarymarker, ...
                                                'InGroup', Inputs.WrapperGroup);

            end
            
            % Seg at bottom
            FemmProblem = addsegments_mfemm(FemmProblem, ...
                                            lastbotnodeid, ...
                                            botnodeid, ...
                                            'BoundaryMarker', FemmProblem.BoundaryProps(boundind).Name, ...
                                            'InGroup', Inputs.WrapperGroup);

            % Add a node at the mid-point of the wrapper, the purpose of this
            % is to ensure an arc segmetn is never more than 180 degrees which
            % causes problems
            [x,y] = pol2cart(thetapole, innerrad);
            [FemmProblem, nodeinds, midnodeid] = addnodes_mfemm(FemmProblem, ...
                                            x, ...
                                            y, ...
                                            'InGroup', Inputs.WrapperGroup);
                                    
            % Seg joining top and bottom (the two most recently added nodes)
            FemmProblem = addarcsegments_mfemm( FemmProblem, ...
                                                botnodeid, ...
                                                midnodeid, ...
                                                rad2deg(thetapole), ...
                                                'InGroup', Inputs.WrapperGroup );

            % Seg joining top and bottom (the two most recently added nodes)
            FemmProblem = addarcsegments_mfemm( FemmProblem, ...
                                                midnodeid, ...
                                                topnodeid, ...
                                                rad2deg(thetapole), ...
                                                'InGroup', Inputs.WrapperGroup );
            
            [innercentres(i,1), innercentres(i,2)] = pol2cart(thetapole, innerrad+wrapperthickness(i,1)/2);
            
        else
            wrapperthickness(i,1) = 0;
        end
        
    end
    
    outercentres = zeros(size(wrapperthickness)) * NaN;
    
    % wrapperthickness(1,2) is the first outer region thickness
    if wrapperthickness(1,2) > Inputs.Tol
        
        outerrad = roffset + rmag/2 + wrapperthickness(1,2);

        % First node is to right of second node in 'nodes' matrix. this is at
        % the bottom of the sim
        [FemmProblem, nodeinds, botnodeid] = addnodes_mfemm(FemmProblem, ...
                                     outerrad, ...
                                     0, ...
                                     'InGroup', Inputs.WrapperGroup);

        if linktb
            lastbotnodeid = elcount.NNodes - size(nodes,1) + 3;
            topnodeid = botnodeid;
            botboundarymarker = '';
        else
            % Second node is to left of penultimate node in 'nodes' matrix
            [x,y] = pol2cart(2*thetapole, outerrad);
            [FemmProblem, nodeinds, topnodeid] = addnodes_mfemm(FemmProblem, ...
                                         x, ...
                                         y, ...
                                        'InGroup', Inputs.WrapperGroup);

            % add a new periodic boundary for the top and bottom of the
            % region
            [FemmProblem, boundind] = addboundaryprop_mfemm(FemmProblem, 'Right Wrap Annular Sec Mags Periodic', 4);

            botboundarymarker = FemmProblem.BoundaryProps(boundind).Name;
            
            % Seg with Periodic boundary at top
            FemmProblem = addsegments_mfemm( FemmProblem, ...
                                             elcount.NNodes - 1, ...
                                             topnodeid, ...
                                             'BoundaryMarker', FemmProblem.BoundaryProps(boundind).Name, ...
                                             'InGroup', Inputs.WrapperGroup );
        
            lastbotnodeid = elcount.NNodes - size(nodes,1) + 1;
        end
        
        % Seg at bottom
        FemmProblem = addsegments_mfemm( FemmProblem, ...
                                         lastbotnodeid, ...
                                         botnodeid, ...
                                         'BoundaryMarker', botboundarymarker, ...
                                         'InGroup', Inputs.WrapperGroup);

        % Add a node at the mid-point of the wrapper, the purpose of this
        % is to ensure an arc segment is never more than 180 degrees which
        % causes problems
        [x,y] = pol2cart(thetapole, outerrad);
        [FemmProblem, nodeinds, midnodeid] = addnodes_mfemm(FemmProblem, ...
                                            x, ...
                                            y, ...
                                            'InGroup', Inputs.WrapperGroup);

        % Seg joining top and bottom (the two most recently added nodes)
        FemmProblem = addarcsegments_mfemm( FemmProblem, ...
                                            botnodeid, ...
                                            midnodeid, ...
                                            rad2deg(thetapole), ...
                                            'InGroup', Inputs.WrapperGroup );

        % Seg joining top and bottom (the two most recently added nodes)
        FemmProblem = addarcsegments_mfemm( FemmProblem, ...
                                            midnodeid, ...
                                            topnodeid, ...
                                            rad2deg(thetapole), ...
                                            'InGroup', Inputs.WrapperGroup );

        [outercentres(1,1), outercentres(1,2)] = pol2cart(thetapole, outerrad-wrapperthickness(1,2)/2);
        
    else
        % Set the region thickness to be exactly zero so this can be tested
        % later
        wrapperthickness(1,2) = 0;
    end
    
    % now add all subsequent right hand wrappers
    for i = 2:size(wrapperthickness, 1)
        
        if wrapperthickness(i,2) > Inputs.Tol
            
            outerrad = outerrad + wrapperthickness(i,2);
            
            lastbotnodeid = botnodeid;
            lasttopnodeid = topnodeid;

            % First node is to right of second node in 'nodes' matrix. this is at
            % the bottom of the sim
            [FemmProblem, nodeinds, botnodeid] = addnodes_mfemm(FemmProblem, ...
                                         outerrad, ...
                                         0, ...
                                         'InGroup', Inputs.WrapperGroup);

            if linktb
                topnodeid = botnodeid;
                botboundarymarker = '';
            else
                % Second node is to left of penultimate node in 'nodes' matrix
                [x,y] = pol2cart(2*thetapole, outerrad);
                [FemmProblem, nodeinds, topnodeid] = addnodes_mfemm(FemmProblem, ...
                                             x, ...
                                             y, ...
                                            'InGroup', Inputs.WrapperGroup);

                % add a new periodic boundary for the top and bottom of the
                % region
                [FemmProblem, boundind] = addboundaryprop_mfemm(FemmProblem, 'Right Wrap Annular Sec Mags Periodic', 4);

                botboundarymarker = FemmProblem.BoundaryProps(boundind).Name;
                
                % Seg with Periodic boundary at top
                FemmProblem = addsegments_mfemm( FemmProblem, ...
                                                 lasttopnodeid, ...
                                                 topnodeid, ...
                                                 'BoundaryMarker', FemmProblem.BoundaryProps(boundind).Name, ...
                                                 'InGroup', Inputs.WrapperGroup );
            end
            
            % Seg at bottom
            FemmProblem = addsegments_mfemm( FemmProblem, ...
                                             lastbotnodeid, ...
                                             botnodeid, ...
                                             'BoundaryMarker', botboundarymarker, ...
                                             'InGroup', Inputs.WrapperGroup );

            % Add a node at the mid-point of the wrapper, the purpose of this
            % is to ensure an arc segmetn is never more than 180 degrees which
            % causes problems
            [x,y] = pol2cart(thetapole, outerrad);
            [FemmProblem, nodeinds, midnodeid] = addnodes_mfemm(FemmProblem, ...
                                            x, ...
                                            y, ...
                                            'InGroup', Inputs.WrapperGroup);

            % Seg joining top and bottom (the two most recently added nodes)
            FemmProblem = addarcsegments_mfemm( FemmProblem, ...
                                                botnodeid, ...
                                                midnodeid, ...
                                                rad2deg(thetapole), ...
                                                'InGroup', Inputs.WrapperGroup );

            % Seg joining top and bottom (the two most recently added nodes)
            FemmProblem = addarcsegments_mfemm( FemmProblem, ...
                                                midnodeid, ...
                                                topnodeid, ...
                                                rad2deg(thetapole), ...
                                                'InGroup', Inputs.WrapperGroup );
            
            [outercentres(i,1), outercentres(i,2)] = pol2cart(thetapole, outerrad-wrapperthickness(i,2)/2);
            
        else
            wrapperthickness(i,2) = 0;
        end
        
    end

end
