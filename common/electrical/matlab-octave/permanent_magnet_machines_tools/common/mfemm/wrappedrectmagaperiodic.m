function [FemmProblem, wrapperthickness, leftcentres, rightcentres, outernodeids] = wrappedrectmagaperiodic(FemmProblem, ypole, ymag, xmag, xoffset, pos, wrapperthickness, varargin)
% draws a two-pole magnets set with 'wrappers' on either side of varying
% thickness and aperiodic boudaries at the top and bottom
%
% Syntax
%
% [FemmProblem, wrapperthickness, leftcentres, rightcentres, outernodeids] = ...
%       wrappedrectmagaperiodic(FemmProblem, ypole, ymag, xmag, xoffset, pos, wrapperthickness, 'Parameter', Value)
%
% Input
%
% 

    Inputs.MagDirections = {0, 180};
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
    
    % remove input fields specific to this function and pass the
    % remaining fields to planarrectmagaperiodic as a series of p-v pairs
    planarrectmagaperiodicInputs = rmfield(Inputs, {'WrapperGroup'});
    
    % Convert the inputs structure to a series of p-v pairs
    planarrectmagaperiodicInputs = struct2pvpairs(planarrectmagaperiodicInputs);
    
    % first draw periodic magnet regions
    [FemmProblem, nodes, nodeids, links] = rectmagaperiodic(FemmProblem, ypole, ymag, xmag, xoffset, pos, planarrectmagaperiodicInputs{:});
    
    leftouternodeids = [nodeids(1), nodeids(end-1)];
    rightouternodeids = [nodeids(2), nodeids(end)];
    
    % now add the back iron nodes and links, depending on their thicknesses
    if size(wrapperthickness,2) < 2
        wrapperthickness = [wrapperthickness, wrapperthickness];
    elseif size(wrapperthickness,2) > 2
        error('wrapperthickness must be a scaler or a (1 x 2) vector or (n x 2) matrix')
    end
    
    if anyalldims(wrapperthickness < 0)
        error('wrapper thicknesses must all be greater that 0')
    end
    
    elcount = elementcount_mfemm(FemmProblem);
    
    leftcentres = zeros(size(wrapperthickness)) * NaN;
    
    % wrapperthickness(1,1) is the first left hand region thickness
    if wrapperthickness(1,1) > Inputs.Tol
        % add the nodes and segments for the left hand side
        
        % First node is to left of first node in 'nodes' matrix. this is at
        % the bottom of the sim
        FemmProblem = addnodes_mfemm(FemmProblem, nodes(1,1) - wrapperthickness(1,1), nodes(1,2), ...
                            'InGroup', Inputs.WrapperGroup);
        
        % Second node is to left of penultimate node in 'nodes' matrix
        FemmProblem = addnodes_mfemm(FemmProblem, nodes(end-1,1) - wrapperthickness(1,1), nodes(end-1,2), ...
                            'InGroup', Inputs.WrapperGroup);
        
        % add a new periodic boundary for the top and bottom of the region
        [FemmProblem, boundind] = addboundaryprop_mfemm(FemmProblem, 'Left Wrap Rect Mags Periodic', 4);
        
        % Seg with Periodic boundary at bottom
        FemmProblem = addsegments_mfemm(FemmProblem, elcount.NNodes - size(nodes,1), elcount.NNodes, ...
            'BoundaryMarker', FemmProblem.BoundaryProps(boundind).Name, ...
            'InGroup', Inputs.WrapperGroup);

        % Seg with Periodic boundary at top
        FemmProblem = addsegments_mfemm(FemmProblem, elcount.NNodes - 2, elcount.NNodes + 1, ...
            'BoundaryMarker', FemmProblem.BoundaryProps(boundind).Name, ...
            'InGroup', Inputs.WrapperGroup);
        
        % Seg joining top and bottom (the two most recently added nodes)
        FemmProblem = addsegments_mfemm(FemmProblem, numel(FemmProblem.Nodes) - 2, numel(FemmProblem.Nodes) - 1, ...
                        'InGroup', Inputs.WrapperGroup);
        
        leftouternodeids = [numel(FemmProblem.Nodes) - 2, numel(FemmProblem.Nodes) - 1];
        
        leftcentres(1,:) = rectcentre(nodes(1,1) - wrapperthickness(1,1), nodes(1,2), ...
                                      nodes(end-1,1), nodes(end-1,2));
        

    else
        % Set the region thickness to be exactly zero so this can be tested
        % later
        wrapperthickness(1,1) = 0;
    end
    
    % now add all subsequent left hand wrappers
    for i = 2:size(wrapperthickness, 1)
        
        if wrapperthickness(i,1) > Inputs.Tol
        
            % First node is to left of second last node added matrix. this is at
            % the bottom of the sim
            FemmProblem = addnodes_mfemm(FemmProblem, ...
                                         FemmProblem.Nodes(end-1).Coords(1) - wrapperthickness(i,1), ...
                                         FemmProblem.Nodes(end-1).Coords(2), ...
                                         'InGroup', Inputs.WrapperGroup);

            % Second node is to left of the last node added
            FemmProblem = addnodes_mfemm(FemmProblem, ...
                                         FemmProblem.Nodes(end-1).Coords(1) - wrapperthickness(i,1), ...
                                         FemmProblem.Nodes(end-1).Coords(2), ...
                                         'InGroup', Inputs.WrapperGroup);

            % add a new periodic boundary for the top and bottom of the region
            [FemmProblem, boundind] = addboundaryprop_mfemm(FemmProblem, 'Left Wrap Rect Mags Periodic', 4);

            % Seg with Periodic boundary at bottom
            FemmProblem.Segments(end + 1) = newsegment_mfemm( numel(FemmProblem.Nodes) - 2, numel(FemmProblem.Nodes) - 4, ...
                'BoundaryMarker', FemmProblem.BoundaryProps(boundind).Name, ...
                'InGroup', Inputs.WrapperGroup);

            % Seg with Periodic boundary at top
            FemmProblem.Segments(end + 1) = newsegment_mfemm(numel(FemmProblem.Nodes) - 1, numel(FemmProblem.Nodes) - 3, ...
                'BoundaryMarker', FemmProblem.BoundaryProps(boundind).Name, ...
                'InGroup', Inputs.WrapperGroup);
            
            % Seg joining top and bottom (the two most recently added nodes)
            FemmProblem.Segments(end + 1) = newsegment_mfemm(numel(FemmProblem.Nodes) - 2, numel(FemmProblem.Nodes) - 1, ...
                'InGroup', Inputs.WrapperGroup);
            
            leftouternodeids = [numel(FemmProblem.Nodes) - 2, numel(FemmProblem.Nodes) - 1];
            
            leftcentres(i,:) = rectcentre(FemmProblem.Nodes(end-1).Coords(1), FemmProblem.Nodes(end-1).Coords(2), ...
                                          FemmProblem.Nodes(end-2).Coords(1), FemmProblem.Nodes(end-2).Coords(2));
        
        else
            wrapperthickness(i,1) = 0;
        end
        
    end
    
    rightcentres = zeros(size(wrapperthickness)) * NaN;
    
    % wrapperthickness(1) is the left hand region thickness
    if wrapperthickness(1,2) > Inputs.Tol
        
        % First node is to right of second node in 'nodes' matrix. this is at
        % the bottom of the sim
        FemmProblem = addnodes_mfemm(FemmProblem, ...
                                     nodes(2,1) + wrapperthickness(1,2), ...
                                     nodes(2,2), ...
                                     'InGroup', Inputs.WrapperGroup);
        
        % Second node is to right of last node in 'nodes' matrix
        FemmProblem = addnodes_mfemm(FemmProblem, ...
                                     nodes(end,1) + wrapperthickness(1,2), ...
                                     nodes(end,2), ...
                                     'InGroup', Inputs.WrapperGroup);
        
        % add a new periodic boundary for the top and bottom of the
        % region
        [FemmProblem, boundind] = addboundaryprop_mfemm(FemmProblem, 'Right Wrap Rect Mags Periodic', 4);
        
        % Seg with Periodic boundary at bottom
        FemmProblem.Segments(end + 1) = newsegment_mfemm( elcount.NNodes - size(nodes,1) + 1, numel(FemmProblem.Nodes) - 2, ...
            'BoundaryMarker', FemmProblem.BoundaryProps(boundind).Name, ...
            'InGroup', Inputs.WrapperGroup);

        % Seg with Periodic boundary at top
        FemmProblem.Segments(end + 1) = newsegment_mfemm( elcount.NNodes - 1, numel(FemmProblem.Nodes) - 1, ...
            'BoundaryMarker', FemmProblem.BoundaryProps(boundind).Name, ...
            'InGroup', Inputs.WrapperGroup);
        
        % Seg joining top and bottom (the two most recently added nodes)
        FemmProblem.Segments(end + 1) = newsegment_mfemm(numel(FemmProblem.Nodes) - 2, numel(FemmProblem.Nodes) - 1, ...
            'InGroup', Inputs.WrapperGroup);
        
        rightouternodeids = [numel(FemmProblem.Nodes) - 2, numel(FemmProblem.Nodes) - 1];
        
        rightcentres(1,:) = rectcentre(nodes(2,1), nodes(2,2), ...
                                       nodes(end,1) + wrapperthickness(1,2), nodes(end,2));
        
    else
        % Set the region thickness to be exactly zero so this can be tested
        % later
        wrapperthickness(1,2) = 0;
    end
    
    % now add all subsequent right hand wrappers
    for i = 2:size(wrapperthickness, 1)
        
        if wrapperthickness(i,2) > Inputs.Tol
        
            % First node is to left of second last node added. this is at
            % the bottom of the sim
            FemmProblem = addnodes_mfemm(FemmProblem, ...
                                         FemmProblem.Nodes(end-1).Coords(1) + wrapperthickness(i,2), ...
                                         FemmProblem.Nodes(end-1).Coords(2), ...
                                         'InGroup', Inputs.WrapperGroup);

            % Second node is to left of the last node added
            FemmProblem = addnodes_mfemm(FemmProblem, ...
                                         FemmProblem.Nodes(end-1).Coords(1) + wrapperthickness(i,2), ...
                                         FemmProblem.Nodes(end-1).Coords(2), ...
                                         'InGroup', Inputs.WrapperGroup);

            % add a new periodic boundary for the top and bottom of the
            % region
            [FemmProblem, boundind] = addboundaryprop_mfemm(FemmProblem, 'Right Wrap Rect Mags Periodic', 4);

            % Seg with Periodic boundary at bottom
            FemmProblem.Segments(end + 1) = newsegment_mfemm( numel(FemmProblem.Nodes) - 2, numel(FemmProblem.Nodes) - 4, ...
                'BoundaryMarker', FemmProblem.BoundaryProps(boundind).Name, ...
                'InGroup', Inputs.WrapperGroup);

            % Seg with Periodic boundary at top
            FemmProblem.Segments(end + 1) = newsegment_mfemm(numel(FemmProblem.Nodes) - 1, numel(FemmProblem.Nodes) - 3, ...
                'BoundaryMarker', FemmProblem.BoundaryProps(boundind).Name, ...
                'InGroup', Inputs.WrapperGroup);

            % Seg joining top and bottom (the two most recently added nodes)
            FemmProblem.Segments(end + 1) = newsegment_mfemm(numel(FemmProblem.Nodes) - 2, numel(FemmProblem.Nodes) - 1, ...
                'InGroup', Inputs.WrapperGroup);
            
            rightouternodeids = [numel(FemmProblem.Nodes) - 2, numel(FemmProblem.Nodes) - 1];
            
            rightcentres(i,:) = rectcentre(FemmProblem.Nodes(end-1).Coords(1), FemmProblem.Nodes(end-1).Coords(2), ...
                                           FemmProblem.Nodes(end-2).Coords(1), FemmProblem.Nodes(end-2).Coords(2));
        
        else
            wrapperthickness(i,2) = 0;
        end
        
    end
    
    outernodeids = [ leftouternodeids(1), rightouternodeids(1), rightouternodeids(2), leftouternodeids(2) ];

end